from optparse import OptionParser
from scipy import interpolate 
from scipy import optimize
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import cPickle as pickle
import os

#this is the usage message for the option parser
usage="""
Reads a file containing output from Aaron's Kerr script and produces
Omega, OmegaMag
"""

def error(msg):
    os.sys.stderr.write("#### ERROR ####\n");
    os.sys.stderr.write(msg)
    os.sys.stderr.write('\n')
    os.sys.exit(-1)

########################################
#command line arguments begin
p=OptionParser(usage=usage)
p.add_option("--input",type="string",
       default="KerrAnalyticOrbitData.dat", 
           help="The input file.")
p.add_option("--output",type="string",
       default="Omega", 
        help="Name of output. A [out].pkl containing the data and [out].pdf of"+
            " a T/OmegaMag plot will be saved. Open the .pkl with "+
            "e.g. Omega,OmegaMag = pickle.load([out])")
p.add_option("--n",type="int", default=10,
        help="Frequencies will be computed from fits based on 2n+1 data points"+
            " around each point in the input.")
p.add_option("--no_save",action="store_true",
        default=False, 
        help="If specified, don't save the output.")
p.add_option("--no_show",action="store_true",
        default=False, 
        help="If specified, don't display a t vs. OmegaMag plot.")
(opts,args) = p.parse_args()

#main
#open the file
############################################
def Compute_OrbitalFrequency(xA, N):
    """Given numpy arrays for the black hole locations xA, xB,
    perform fits to 2N+1 data-points around each point, and from
    the fit compute the instantaneous orbital frequency,
            Omega = r\times \dot{r} / r^2
    return Omega, OmegaMag
"""
    # compute orbital frequency as a function of time
    fitfunc = lambda p,x: p[0]+p[1]*x+p[2]*x**2
    errfunc = lambda p,x,y: fitfunc(p,x)-y

    def FitData(data, N):
        """given a numpy array data with time as first column, perform fits 
        covering N points before and after each data-point for each column.
        return the fitted values and their first time-derivatives as a 
        numpy arrays, with first column time"""
        # collect output data, first column being time
        last_idx=len(data)-N-1
        x_fit=[ list(data[N:last_idx,0]) ]
        v_fit=[ list(data[N:last_idx,0]) ]
        for ycol in range(1,data.shape[1]):
            p0=[data[N,ycol],0,0]
            x_tmp=[]
            v_tmp=[]
            for idx in range(N, len(data)-N-1):
                T=data[idx,0]
                x=data[idx-N:idx+N+1,0]-T # shift target to zero
                y=data[idx-N:idx+N+1,ycol]
                p0,success=optimize.leastsq(errfunc, p0[:],args=(x,y))
                x_tmp.append(p0[0])
                v_tmp.append(p0[1])
            x_fit.append(x_tmp)
            v_fit.append(v_tmp)
        return np.transpose(np.array(x_fit)),\
               np.transpose(np.array(v_fit))

    xA_fit,vA=FitData(xA,N)
    #compute Orbital frequency
    dr=xA_fit
    dv=vA
    r2=(dr[:,1]**2+dr[:,2]**2+dr[:,3]**2)
    Omega = xA_fit
    Omega[:,1]=(dr[:,2]*dv[:,3]-dr[:,3]*dv[:,2])/r2
    Omega[:,2]=(dr[:,3]*dv[:,1]-dr[:,1]*dv[:,3])/r2
    Omega[:,3]=(dr[:,1]*dv[:,2]-dr[:,2]*dv[:,1])/r2
    OmegaMag=Omega[:,0:2].copy()
    OmegaMag[:,1]=np.sqrt(Omega[:,1]**2+Omega[:,2]**2+Omega[:,3]**2)
    return Omega, OmegaMag
###################################################################

with open(opts.input,'r') as myfile:
  xA = np.loadtxt(myfile)
Omega,OmegaMag = Compute_OrbitalFrequency(xA,opts.n)
r2 = xA[:,1]**2+xA[:,2]**2+xA[:,3]**2
r = np.sqrt(r2)
t = xA[:,0]
plt.plot(t,r)
plt.show()
print Omega
#output
if(opts.no_save==False):
    print "Saving..."
    pklname = opts.output+".pkl"
    g = open(pklname,'w')
    pickle.dump([Omega,OmegaMag],g)
    g.close()
    print "Done."

if(opts.no_show==False):
    plt.plot(OmegaMag[:,0],OmegaMag[:,1])
    plt.xlabel("Time")
    plt.ylabel("Frequency")
    pltname = opts.output+".pdf"
    plt.savefig(pltname)
    plt.show()
    

