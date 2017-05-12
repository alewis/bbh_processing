"""
Functions to collect data from many directories into a central file.
"""


import numpy as np
from ecctools import make_ecc_dict
from utilities import turning_points
import matplotlib.pyplot as plt
import cPickle as pickle
from os.path import expanduser, isfile




#Convert the ecc_dict into a numpy array
def array_from_ecc_dict(ecc_dict):
  return np.asarray([(k[0],k[1],v) for k, v in ecc_dict.iteritems()])


#Write a header to the central file
def construct_header(headerstrings):
  the_header="" 
  i=0
  for hstr in headerstrings:
    the_header+="# ["+str(i+1)+"] = "+ hstr+"\n" 
    i = i+1
  return the_header


def write_to_central_file(
    datapaths = "~/.bbh/.datapaths",
    eccfilepath  = "~/.bbh/.eccfilepath",
    headerstrings = "~/.bbh/.headerstrings",
    sepfile = "HorizonSepMeasures.dat",
    pname = "Omega.pkl",
    mode = "coord"):
    """
    Given a list of paths to data 'datapaths', this runs a function 'func' 
    on each path, which should return a list of tuples.
    """

def write_central_file(
    datapaths = "~/.bbh/.datapaths",
    eccfilepath  = "~/.bbh/.eccfilepath",
    headerstrings = "~/.bbh/.headerstrings",
    sepfile = "HorizonSepMeasures.dat",
    pname = "Omega.pkl",
    mode = "coord"):


  #Turn tildes into path to home directory
  expanded_path=expanduser(eccfilepath).strip()
  #The path to the central file should be the first line in 
  #central_file_location_path
  outfilename=expanduser(open(expanded_path,'r').read()).strip()
  outpicklename = outfilename+".pkl"
  outdatname = outfilename+".dat"
  #Read in the old dictionary, if it is there
  #This maps r_a/Omega_a tuples to eccentricities. 
  ecc_dict = dict() 
  if isfile(outpicklename):
    try:
      ecc_dict = pickle.load(open(outpicklename,"r"))
    except EOFError:
      print "dict file", outpicklename, " exists but was empty."
  
  #Read in the data paths 
  expanded_data_path = expanduser(datapaths).strip()
  with open(expanded_data_path,'r') as f:
    datapathnames = f.read().splitlines()
  for path in datapathnames:
    path = expanduser(path).strip()
  
  for datapathname in datapathnames:
    print datapathname 
    sepfilename = datapathname+"/"+sepfile
    pfilename = datapathname+"/"+pname
    
    new_dict=make_ecc_dict(sepfilename,pfilename,mode)
    
    #add to the old dictonary; overlapping keys are taken from the new data
    ecc_dict.update(new_dict)

  #Get the header string
  expandheader = expanduser(headerstrings).strip() 
  hstring = open(expandheader,"r").read().splitlines()
  header = construct_header(hstring)
 
  picklefile = open(outpicklename,"w")
  pickle.dump(ecc_dict,picklefile)
  picklefile.close() 

  ecc_array = array_from_ecc_dict(ecc_dict)
  outfile = open(outdatname,"w")
  outfile.write(header)
  np.savetxt(outfile,ecc_array)

write_central_file()
