import os
import numpy 
import scipy.optimize
from optparse import OptionParser
import cPickle as pickle

usage="""
Reads a dat file containing
   t OmegaMag
and attempts to extract periastron advance K for the range of W parameters specified.  Generates an output file W*.dat for each value of W.
"""
 

def error(msg):
    os.sys.stderr.write("#### ERROR ####\n");
    os.sys.stderr.write(msg)
    os.sys.stderr.write('\n')
    os.sys.exit(-1)

################################################################

def FittingFunc(p,
                t,
                OmegaMag,
                T):
    if(p.shape[0]==3):
        return p[0]*(p[1]-t)**p[2] \
            - OmegaMag
    if(p.shape[0]==6):
        return p[0]*(p[1]-t)**p[2] \ + p[3]*numpy.cos(p[4]+p[5]*(t-T) ) \
            - OmegaMag

    return p[0]*(p[1]-t)**p[2] \
        + p[3]*numpy.cos(p[4]+p[5]*(t-T)+p[6]*(t-T)**2) \
        - OmegaMag

################################################################
                

p=OptionParser(usage=usage)
p.add_option("--Omega", type="string",
             help="File containing Omega vs. t (orbital frequency)")
p.add_option("--cos1", action="store_true",
             default=False,
             help="If specified, use a linear function inside the argument.  If not specified, use a quadratic function")
p.add_option("--tmin", type="float",
             help="If specified, do not use any data for t<tmin.  This restricts for which Omega the periastron advance is computed")
p.add_option("--bad_fit_cutoff", type="float",
              default=0.3,
              help="Compute ratio [rms of fit]/[Amplitude of Oscillatory part]. If this ratio is O(1) or larger, than this indicates a bad fit.  Prevent output if the ratio is larger than the given number.  DEFAULT = 0.3")
p.add_option("--verbose", action="store_true",
             default=False,
             help="More verbose output")
p.add_option("--mass",type="float",
             default=1,
             help="Total mass of binary system")
p.add_option("--min",type="float",
             help="Minimum W parameter")
p.add_option("--max",type="float",
             help="Maximum W parameter")
p.add_option("--dw",type="float",
             help="How much to increment W each loop iteration")


(opts,args)=p.parse_args()
if(opts.Omega == None):
    error("Must specify --Omega")
if opts.min is None or opts.max is None or opts.dw is None:
    error("Must specify the limits and increment for W (--min, --max, --dw)")
W=opts.min

while W<(opts.max+opts.dw):

    print 'W=', W
    f = open(opts.Omega,'r') 
    bla,blamag=pickle.load(f)
    t=bla[:,0]
    Omega=blamag[:,1]
    tmin=t[0]
    if(opts.tmin!=None):
        tmin=max(opts.tmin, tmin)

    firstFit=True
    #print t.shape, Omega.shape
    #print t

    # initialize initial parameters of fit 
    if(opts.cos1):
        p =numpy.arange(0.,6)
    else:
        p =numpy.arange(0.,7)
        p[6]=0.;

    p[0]=0.4
    p[1]=t[-1]+1500.
    p[2]=-0.35
    p[3]=0.  # oscillatory part initialized during firstFit below
    p[4]=0.
    p[5]=0.
    Params=[]
    Ts=[]
    RMSs=[]
    for idx in range(0,t.shape[0],10):
        T=t[idx]
        O=Omega[idx]
        if O<0:
          raise ValueError("Encountered negative frequency; aborting")

        
        width= W*2*numpy.pi / O
        if(T-width < tmin) or (T+width > t[-1]):
            print "Skipping T=%g (lack of data)"%T
            continue
        cut_indices=numpy.abs(t-T)<width
        tcut=t[ cut_indices ]
        Ocut=Omega[ cut_indices ]
        #print tcut.shape
        #print tcut
        if(firstFit):
            pin=p[0:3]
            if(opts.verbose):
                print "First fit, start with preparatory F1cos0"
                print "   - pin= ", pin
            pout, cov_v, infodict, msg, ier=scipy.optimize.leastsq(FittingFunc, pin, args=(tcut, Ocut, T), full_output=1)

            residual=FittingFunc(pout, tcut,Ocut,T)
            rms=numpy.sqrt(numpy.mean(residual**2))

            # now construct full initial guess for leastsq
            #  - smooth part from preparatory fit
            #  - amplitude of oscillatory part from rms of preparatory fit
            #  - phase based on residual in middle at T
            #  - Omega_r based on Omega_phi
            p[0:3]=pout
            p[3]=rms* 2**0.5  # sqrt(2) converts from rms to peak
            if(numpy.mean(residual[numpy.abs(tcut-T)<5.]) < 0):
                p[4]=0.
            else:
                p[4]=numpy.pi 
            p[5] = O/1.3
            if(opts.verbose):
                print "   - pout=",pout
                print "   - pnew=",p
            firstFit=False

        #print "pin= ",p
        pout, cov_v, infodict, msg, ierr=scipy.optimize.leastsq(FittingFunc, p, args=(tcut, Ocut, T), full_output=1)
        rms=numpy.sqrt(numpy.mean(FittingFunc(p,tcut,Ocut,T)**2))
        if(opts.verbose):
            print "T=%6.1f, W=%6.1f: ierr=%i, nfunc=%i; rms=%4e, rms/B=%g"%(T, width, ierr, infodict['nfev'], rms, rms/numpy.abs(pout[3]))
        #print "Output  parameters for fit pout=",pout
        # use this optimization as new guess
        p=pout

        # remember for output
        if(rms/numpy.abs(pout[3]) > opts.bad_fit_cutoff):
            print "SKIPPING OUTPUT at T=%g.  Reason: residual/B > bad_fit_cutoff."%T
            continue
        if(pout[5]<0):
            print "SKIPPING OUTPUT at T=%g:  Reason: Fitted Omega_r=%g negative."%(T,pout[5])
            continue


        RMSs.append(rms)
        Params.append( pout)
        Ts.append(T)


    Params=numpy.array(Params)
    Ts=numpy.array(Ts)
    RMSs=numpy.array(RMSs) # root-mean-square residual of fit
    p0=Params[:,0]
    p1=Params[:,1]
    p2=Params[:,2]
    p3=Params[:,3]
    p4=Params[:,4]
    p5=Params[:,5]
    p6=Params[:,6]  #don't output p6, because it isn't always present
    OmegaOrbit=p0*(p1-Ts)**p2
    Omegar=p5
    K=OmegaOrbit/Omegar


    # Schwarzschild
    #    K_S = 1 + 3/r
    # r by Keplers law
    #    r^3 Omega^2 = M,   M/r = (M Omega)^(2/3) 
    # therefore, 
    #    K_S = 1 + 3 (M Omega)^(2/3)
    K_over_KS = K *numpy.sqrt (1-6*(OmegaOrbit*opts.mass)**(2./3))

    e_OmegaPhi = numpy.abs(p3)/(2.*OmegaOrbit)
    e_Omegar=numpy.abs(p3)/(2.*Omegar)

    RMS_over_B = RMSs/p3

    out=file("W"+str(W)+".dat",'w')
    out.write("""#
    # [1] = t
    # [2] = Omega_phi
    # [3] = Omega_r
    # [4] = K
    # [5] = K/K_S
    # [6] = abs(p3)/2Omega_phi
    # [7] = abs(p3)/2Omega_r
    # [8] = residual of fit
    # [9] = residual/B
    # [10] = p0
    # [11] = p1
    # [12] = p2
    # [13] = p3
    # [14] = p4
    # [15] = p5
    # [16] = p6
    """)

    numpy.savetxt(out,
                  numpy.transpose([Ts,opts.mass*OmegaOrbit,opts.mass*Omegar,K,K_over_KS, 
                                   e_OmegaPhi, e_Omegar, RMSs, RMS_over_B, p0,p1,p2,p3,p4,p5,p6]))

    W+=opts.dw
