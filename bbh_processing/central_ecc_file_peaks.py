import numpy as np
from ecctools import make_ecc_peak_array 
from utilities import turning_points
import matplotlib.pyplot as plt
import cPickle as pickle
from os.path import expanduser, isfile

#Save a .pkl file storing a list of numpy arrays,
#each one storing separation,frequency,and ecc data
#from each simulation

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
  
  #Read in the data paths 
  expanded_data_path = expanduser(datapaths).strip()
  with open(expanded_data_path,'r') as f:
    datapathnames = f.read().splitlines()
  for path in datapathnames:
    path = expanduser(path).strip()

  arraylist=[]
  for datapathname in datapathnames:
    print datapathname 
    sepfilename = datapathname+"/"+sepfile
    pfilename = datapathname+"/"+pname
    try:
      errcode, this_array = make_ecc_peak_array(sepfilename,pfilename,mode,"interp") 
      if (errcode==0): 
        arraylist.append(this_array)
      else:
        print "Received error code: Not enough peaks for this run to fit spline"
    except IOError:
      print "IOError: file missing for this path."

  picklefile = open(outpicklename,"w")
  pickle.dump(arraylist,picklefile)
  picklefile.close() 
  

write_central_file()
