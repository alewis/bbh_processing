"""Contains methods to interact with Aaron's Kerr integrator and to manipulate
the results.
"""
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from os.path import isfile
import analyze_kerr_data as kerr
import string

def generate_kerr_orbit(r_a, omega_0, inc, spin, n_orbits, script_path,
                        out_path):
    """Run Aaron's Kerr integrator to generate a bound Kerr orbit with the
    given parameters. Params:
          r_a -> apostron radius (usually D0)
          omega_0 -> initial orbital frequency
          inc -> maximal orbital inclination in degrees
          spin -> Kerr spin parameter
          n_orbits -> the number of orbits to integrate over
          script_path -> absolute path to WriteKerrOrbit.m
          out_path -> path to location of output
    You should have Mathematica >=10 installed with its module loaded (so that
    entering 'WolframScript' at the command line has an effect). 
    """
    print "Input params were: ",
    print "r_a = ", r_a, "omega_0 = ", omega_0, "inc = ", inc 
    print "spin = ", spin, "n_orbits = ", n_orbits
    print "script_path =", script_path
    
    tag_string = "KerrOrbit_rA_"+str(r_a)+"_om_"+str(omega_0)
    tag_string += "_inc_"+str(inc)+"_a_"+str(spin)+"_norb_"+str(n_orbits)
    tag_string = tag_string.replace(".", "_") 
    out_string = out_path + tag_string
    if not isfile(script_path):
        raise ValueError("ERROR: no file found at location\n "+script_path)
    callstring = "WolframScript -script "+script_path
    callstring += " "+str(r_a)+" "+str(omega_0)+" "+str(inc)+" "
    callstring += str(spin) + " " + str(n_orbits) + " " + str(out_string)
    print "Shell call: ", callstring
    call([callstring], shell=True)
    print "Call succeeded!"
    return out_string + ".dat"

if __name__=="__main__":
    pathstring = "/mnt/raid-project/nr/adlewis/ResonanceAnalysis/Notes/"
    pathstring += "GenericBBH/IPythonNotebooks/WriteKerrOrbit.m"
    outstring = "./orbits/"
    output = generate_kerr_orbit(10, 0.004, 20, 0.8, 2, pathstring, outstring)
    print output
    kerr.loadfile(output)
