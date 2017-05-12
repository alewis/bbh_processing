import numpy as np
from subprocess import Popen 
from os.path import isfile
import pandas as pd
import itertools
def tostr(n):
    return str(n).replace(".", "_")

def make_filename(p, ecc, zm, a, norbits, path):
    fname = "/KerrOrbit_p"+tostr(p)+"_Ecc"+tostr(ecc)+"_Zm"+tostr(zm)
    fname += "_Spin"+tostr(a)+"_Norbits"+tostr(norbits)+".dat"
    return path+fname

def make_callstring(p, ecc, zm, a, norbits, path, script_path):
    outpath = make_filename(p, ecc, zm, a, norbits, path)
    print outpath
    callstring = "WolframScript -script "+script_path
    callstring += " "+str(p)+" "+str(ecc)+" "+str(zm)+" "
    callstring += str(a) + " " + str(norbits) + " " + outpath 
    return callstring

def generate_files(callstrs, script_path):
    """Run Aaron's Kerr integrator to generate a bound Kerr orbit with the
    given parameters. :
    You should have Mathematica >=10 installed with its module loaded (so that
    entering 'WolframScript' at the command line has an effect). 
    """
    if not isfile(script_path):
        raise ValueError("ERROR: no file found at location\n "+script_path)
    for callstring in callstrs:
        print "Shell call: ", callstring
        p1 = Popen([callstring], shell=True)
        p1.wait()
    print "Call succeeded!"

def kerr_dataframes(plist, ecclist, zmlist, spinlist, norblist, outpath=None,
                    scriptpath = None):
    """ Make a dataframe of kerr orbits. These have semi-latus recta, 
        eccentricities, inclinations, and spins given by the arguments above.
        Outpath points to the location we will be saving the orbits to.
        Scriptpath points to the location of the Wolfram script running Aaron's
        
    """
    if outpath==None:
      outpath = "/cita/d/raid-project/nr/adlewis/ResonanceAnalysis/Notes/"
      outpath += "GenericBBH/ExactOrbits/orbits"
    
    joined = itertools.product(plist, ecclist, zmlist, spinlist, norblist)
    fnames = [make_filename(p, ecc, zm, a, norb, outpath) 
              for p, ecc, zm, a, norb in joined]
    joined = itertools.product(plist, ecclist, zmlist, spinlist, norblist)
    callstrings = [make_callstring(p, ecc, zm, a, norb, outpath, scriptpath) 
              for p, ecc, zm, a, norb in joined]
    generate_files(callstrings, scriptpath)
    entries = [read_orbit(fname) for fname in fnames] 
    the_df = pd.DataFrame(entries)
    return the_df

def read_orbit(fname):
    """ Reads one of Aaron's orbit files and returns the data as a Series.
    """
    traj = np.loadtxt(fname)
   
    #Read the header
    the_dict = {"traj": traj}
    with open(fname, "r") as f:
        for line in f:
            if line[0] != "#":
                break
            else:
                split = line[1:].split()
                tag = split[0]
                dat = float(split[-1])
                the_dict[tag] = dat
    the_series = pd.Series(the_dict)
    return the_series

def make_orbits():
    eccmin = 0
    eccmax = 0.9
    r_pmin = 6.
    r_pmax = 15
    
    eccsqlist = np.linspace(1-eccmin**2, 1-eccmax**2, 10)
    ecclist = np.sqrt(1-eccsqlist)
    rplist = np.linspace(r_pmin, r_pmax, 10) 
    
    rho = lambda rp, ecc: rp*(1-ecc**2)/(1-ecc)
    plist = rho(rplist, ecclist) 
    
    zmlist = np.linspace(0, 0.9, 10)
    spinlist = [0.8]
    norblist = [10]

    stem = "/cita/d/raid-project/nr/adlewis/ResonanceAnalysis/Notes/GenericBBH"
    outpath = stem+"/ExactOrbits/orbits"

    scriptpath = stem+"/ExactOrbits/WriteKerrOrbit.m"
    df = kerr_dataframes(plist, ecclist, zmlist, spinlist, norblist, outpath,
                    scriptpath)
    df.to_pickle(outpath+"/orbit_dfs.pkl")          

if __name__ == "__main__":
    make_orbits()    

