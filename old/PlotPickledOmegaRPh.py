import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
from optparse import OptionParser

#this is the usage message for the option parser
usage="""
Reads output from GenerateOmegaRPh.py, plots it, and saves and image of the 
result. 
"""

prs=OptionParser(usage=usage)
prs.add_option("--target",type="string",
        help="Name of the .pkl file to load without the .pkl extension.")
prs.add_option("--magfile",type="string",
        help="Name of the .pkl file containing the total frequencies")
prs.add_option("--outputname", type="string", default="usetarget",
    help="Name of the file to save as, with extension. File format can be "+
        "changed by changing the extension. Default: (target).pdf")
(opts,args) = prs.parse_args()

fname = opts.target+".pkl"
f=open(fname,'r')
T,Omega_ph,Omega_r,totalfit,residuals = pickle.load(f)
f.close()

fname = opts.magfile
f = open(fname,'r')
Omega, OmegaMag = pickle.load(f)
f.close()

plt.figure()
plt.plot(T,Omega_r,label='Om_r')
plt.plot(T,Omega_ph,label='Om_ph')
plt.plot(OmegaMag[:,0],OmegaMag[:,1],label='Fitted data')
plt.plot(T,totalfit,label='Total fit')
plt.legend()
plt.xlabel("Time")
plt.ylabel("Frequency")

if (opts.outputname=="usetarget"):
    plotname=opts.target+".pdf"
else:
    plotname = opts.outputname
plt.title(plotname)
plt.savefig(plotname)
plt.show(plotname)
