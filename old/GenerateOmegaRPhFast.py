import h5py
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
from optparse import OptionParser
import os
#import cProfile as profile
#this is the usage message for the option parser
usage="""
Reads output from GenerateOmegaTotal.py and fits to OmegaMag in such
a way that radial and phi frequencies can be separated. Outputs pickled
    T,Omega_ph,Omega_r,totalfit,residuals
with T time, Omega_i the i frequency, totalfit the full computed fit, and
residuals the (absolute) difference between this fit and OmegaMag. 
"""

def error(msg):
    os.sys.stderr.write("#### ERROR ####\n");
    os.sys.stderr.write(msg)
    os.sys.stderr.write('\n')
    os.sys.exit(-1)
########################################
#command line arguments begin
prs=OptionParser(usage=usage)
prs.add_option("--fit_type",type="string",
        help="Labels the desired fitting function " +
         "(but note this is further modified by the argument --quadratic_term "
         +"defined below).\n"+
         "All cases involve a monotonic term \n"+
         " \t Om(t)_mon = p0* (p1 - (t-T))**p2 \n"+
         "and an oscillatory term \n "+
         " \t Om(t)_osc = p3*cos(p4 + p5*(t-T) + p6*(t-T)**2) \n"+
         " with the p6 term included only if --quadratic_term is specified. \n"+
         " \r--fit_type='add_osc': \n"+
         "\t Om(t) = Om(t)_mon + Om(t)_osc \n"+
         " --fit_type='divide_osc' \n"+
         "\t Om(t) = Om(t)_mon / (1-Om(t)_osc). ")
prs.add_option("--quadratic_term",action="store_true",
        default=False,
        help="If specified, include a (t-T)**2 term in the fit; see above.")
prs.add_option("--input",type="string",
       default="Omega.pkl", 
           help="GenerateOmegaTotal.py output holding Omega and OmegaMag. "+ 
           "DEFAULT:Omega.pkl")
prs.add_option("--output",type="string",
       default="Omega", 
        help="Name of output **without extension**. An [out]+str(n).pkl " +
        " containing the data and an [out]+str(n).pdf plot will be saved. "+
        " Open the .pkl with e.g. \n "+
        "T,Omega_ph,Omega_r,totalfit,residuals=pickle.load([out]+str(n)+'.pkl'). "+ 
        "\n DEFAULT: OmegaRPh")
prs.add_option("--n",type="int", default=2000,
        help="The window of data about each point to sample for fitting. " +
            " Should contain multiple orbits. ")
prs.add_option("--no_save",action="store_true",
        default=False, 
        help="If specified, don't save the output.")
prs.add_option("--no_show",action="store_true",
        default=False, 
        help="If specified, don't save a t vs. OmegaMag plot.")

(opts,args) = prs.parse_args()
############################################


#define possible fitting functions
#avoid nested function calls for speed


#the monotonic piece - always used
om_mon = lambda p,t,T,y: p[0]*((p[1]-(t-T))**p[2]) - y
om_mon_tequalsT = lambda p: p[0]*((p[1])**p[2])
#the oscillatory piece may have a quadratic term
if(opts.quadratic_term==False):
    print "Using linear term inside oscillatory piece."
    om_osc = lambda p,t,T,y: p[0]*np.cos(p[1]+p[2]*(t-T)) - y
else:
    print "Using quadratic term inside oscillatory piece."
    om_osc = lambda p,t,T,y: p[0]*np.cos(p[1]+p[2]*(t-T)+p[3]*((t-T)**2)) - y

#the full fit -> avoid nested function calls for speed
if(opts.fit_type=="add_osc"):
    if(opts.quadratic_term==False):
        def om_fit(p,t,T,y):
            tm = t-T
            mon = p[0]*((p[1]-tm)**p[2])
            poly = p[4]+p[5]*tm
            osc = p[3]*np.cos(poly)
            return mon+osc-y
    else:
        def om_fit(p,t,T,y):
            tm = t-T
            mon = p[0]*((p[1]-tm)**p[2])
            poly = p[4]+tm*(p[5]+p[6]*tm)
            osc = p[3]*np.cos(poly)
            return mon+osc-y

elif(opts.fit_type=="div_osc"):
    if(opts.quadratic_term==False):
        def om_fit(p,t,T,y):
            tm = t-T
            mon = p[0]*((p[1]-tm)**p[2])
            poly = p[4]+p[5]*tm
            osc = p[3]*np.cos(poly)
            return (mon/(1-osc))-y
    else:
        def om_fit(p,t,T,y):
            tm = t-T
            mon = p[0]*((p[1]-tm)**p[2])
            poly = p[4]+tm*(p[5]+p[6]*tm)
            osc = p[3]*np.cos(poly)
            return (mon/(1-osc))-y
    
else:
    error("--fit_type must be specified and either 'add_osc' or 'div_osc'.")


def FitOmegaToModel(OmegaMag, delT):
  
   def FitData(OmegaMag, delT):
        # collect output data, first column being time
        delT2 = delT/2
        last_T=len(OmegaMag)-delT2
        last_T = last_T
        
        
        #prepare the output arrays
        #static numpy array are faster so use them instead of append
        outsize = last_T - delT2
        Tarr = OmegaMag[delT2:last_T,0]
        Omega_ph = np.zeros(outsize) 
        Omega_r = np.zeros(outsize)
        totalfit= np.zeros(outsize)
        residuals= np.zeros(outsize)
        
        
        
        #do early fits to make initial guesses of parameters
        yin = OmegaMag[0:delT,1]
        T = OmegaMag[delT2,0]
        tin = OmegaMag[0:delT,0]
        ymax = max(yin)
        ymin = min(yin)
        amp = (ymax-ymin)/2.0 
        om = 2.*3.14159*0.0017
        
        pmon = [1,OmegaMag[delT,0],1]
        

        if(opts.quadratic_term==False):
            posc = [amp,3,om]
        else:
            posc = [amp,3,om,om]
        
        pmon, success = optimize.leastsq(om_mon,pmon[:],args=(tin,T,yin)) 
        print "Pmon initial: ", pmon
        posc, success = optimize.leastsq(om_osc,posc[:],args=(tin,T,yin))
        print "Posc initial: ", posc
        
        #join the guessed parameters into a single array
        ptot=[]
        for p_i in pmon:
            ptot=np.append(ptot,p_i)
        for p_i in posc:
            ptot=np.append(ptot,p_i) 
        print "steps_remaining,", "T,", "Om_ph,", "Om_r","residual"
        
        #loop over timesteps and do the actual fitting
        #pr = profile.Profile()
        for i in range(delT2, last_T):
             
            #pr.enable()
            ib = i-delT2
            it = i+delT2
            t=Tarr[ib]
            #Tarr[ib] = OmegaMag[i,0] #the current timestep
            
            ptot,success=optimize.leastsq(om_fit, ptot[:],
                    args=(OmegaMag[ib:it,0],t,
                        OmegaMag[ib:it,1]))
            #Tarr[ib]=T
            Omega_ph[ib] = ptot[0]*ptot[1]**ptot[2]
            Omega_r[ib] = ptot[5]
            totalfit[ib] = om_fit(ptot[:],t,t,0)
            residuals[ib] = np.absolute(OmegaMag[i,1]-totalfit[ib])
            
            print last_T - i, t, Omega_ph[ib], Omega_r[ib], residuals[ib]
            #pr.disable()
            #pr.print_stats()
        #loop over timesteps and do the actual fitting
        return Tarr, Omega_ph, Omega_r,totalfit,residuals
   Tarr,Omega_ph,Omega_r,totalfit,residuals = FitData(OmegaMag,delT)
   return Tarr, Omega_ph, Omega_r,totalfit,residuals


#main
#extract the OmegaMag data
f = open(opts.input,'r')
Omega, OmegaMag = pickle.load(f)
f.close()
#do the fit
T,Omega_ph,Omega_r,totalfit,residuals = FitOmegaToModel(OmegaMag,opts.n)

if(opts.no_save==False):
    if(opts.quadratic_term==True):
        pklname="OmegaRPh_quadratic_"+opts.fit_type+"_n="+str(opts.n)+".pkl"
    else:
        pklname="OmegaRPh_linear_"+opts.fit_type+"_n="+str(opts.n)+".pkl"
    g=open(pklname,'w')
    pickle.dump([T,Omega_ph,Omega_r,totalfit,residuals],g)
    g.close()

if(opts.no_show==False):
    plt.figure()
    plt.plot(T[:],Omega_r[:],label='Om_r')
    plt.plot(T[:],Omega_ph[:],label='Om_ph')
    plt.plot(OmegaMag[:,0],OmegaMag[:,1],label='Fitted data')
    plt.plot(T[:],totalfit[:],label='Total fit')
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Frequency")
    if(opts.quadratic_term==True):
        plotname="OmegaRPh_quadratic_"+opts.fit_type+"_n="+str(opts.n)+".pdf"
    else:
        plotname="OmegaRPh_linear_"+opts.fit_type+"_n="+str(opts.n)+".pdf"
    plt.title(plotname)
    plt.savefig(plotname)

print "Done."
