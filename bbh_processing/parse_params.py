from string import whitespace
import math
import utilities as bbh
import hfile_tools as hf
import GenerateOmegaTotal as om
import cPickle as pickle
from os.path import isfile
def parse_params(path="./Params.input", hpath="./Horizons.h5", 
                 omegapath="./Omega.pkl", make_omega=False):
    """Generates a dictionary mapping simulation quantities to their values. 
        First, the Params.input file is parsed and the data added. Then 
        certain data are extracted from the Horizons.h5 file. Finally,
        orbital frequency data are taken from Omega.pkl. If make_omega is True 
        we first generate Omega.pkl.
    """
    with open(path, "r") as f:
        lst=[]
        for line in f:
            li = line.strip()
            if not li.startswith("#"):
                if not li == '':
                    lst.append(li)
    #The non-spin entries 
    keys = [l.split()[0][1:] for l in lst[0:4]] #the param names
    keys += [l.split()[0][1:] for l in lst[6:8]]
    vals = [float(l.split()[2][:-1]) for l in lst[0:4]]
    vals += [float(l.split()[2][:-1]) for l in lst[6:8]]

    #The spin entries
    keys += [l.split()[0][1:] for l in lst[4:6]]
    newvals = [l.partition("(")[2].partition(")")[0] for l in lst[4:6]]
    newvals = [n.translate(None, whitespace).split(",") for n in newvals]
    newvals = [tuple([float(i) for i in n]) for n in newvals]
    vals += newvals

    #The horizons.h5 data
    the_dict = dict(zip(keys,vals))
    the_dict = dict_for_ID_analysis(dict(zip(keys,vals))) 
    the_dict = add_separation(the_dict, hpath)  
    
    if not isfile(omegapath):
        print "No ", omegapath, "found; generating..."
        om.generate_frequency(hpath, omegapath, 10, no_save=False, no_show=True)    
    with open(omegapath) as f:
        Omega, OmegaMag = pickle.load(f)
    the_dict["OmegaMag"] = OmegaMag

    return the_dict 

def getinc(spin):
    """
    Compute the inclination angle in radians of spin with the z-axis.
    """
    #magnitude of projection into XY plane
    xymag = math.sqrt(spin[0]**2 + spin[1]**2)
    angle = math.atan2(xymag, spin[2])
    return math.degrees(angle)

def dict_for_ID_analysis(the_dict):
    """
    Add some extra quantities:
    EccN - Newtonian eccentricity prediction
    ChiA - magnitude of SpinA
    ChiB - magnitude of SpinB
    IncA - angle between SpinA and z axis
    IncB - angle between SpinB and z axis
    Eta  - symm. mass ratio
    """
    pairs = []

    #Eta
    q = the_dict["MassRatio"]
    eta = q/(q+1)**2
    the_dict["Eta"] = eta

    #EccN
    d0 = the_dict["D0"]
    omega0 = the_dict["Omega0"]
    ecc = 1.-omega0**2 * d0**3
    the_dict["EccN"] = ecc

    #ChiA
    spinA = the_dict["SpinA"]
    chiA = 0
    for i in spinA:
        chiA += i**2
    chiA = math.sqrt(chiA)
    the_dict["ChiA"] = chiA
    
    #IncA
    incA = getinc(spinA) 
    the_dict["IncA"] = incA

    #ChiB
    spinB = the_dict["SpinB"]
    chiB = 0
    for i in spinB:
        chiB += i**2
    chiB = math.sqrt(chiB)
    pairs += [("ChiB", chiB)]
    the_dict["ChiB"] = chiB
   
    #Inc B
    incB = getinc(spinB) 
    pairs += [("IncB", incB)]
    the_dict["IncB"] = incB
    return the_dict 

def add_separation(the_dict, hpath="./Horizons.h5"):
    """ 
    Adds t vs. r data to the dictionary, extracted from Horizons.h5".
    """
    t = hf.get_coord("t", hpath)
    rmag = hf.get_coord("r", hpath)
    costheta = hf.get_coord("costheta", hpath, rotate=False)
    thetaphase = hf.get_coord("theta_phase", hpath, rotate=False)

    cosphi = hf.get_coord("cosphi", hpath, rotate=True)
    phiphase = hf.get_coord("phi_phase", hpath, rotate=True)

    #rotate into correct plane
    chi = hf.chi_inertial(hpath)
    the_dict["Time"] = t 
    the_dict["r_mag"] = rmag
    the_dict["CosTheta"] = costheta
    the_dict["ThetaPhase"] = thetaphase
    the_dict["CosPhi"] = cosphi
    the_dict["PhiPhase"] = phiphase
    return the_dict 


if __name__ == "__main__":
    parse_params()
