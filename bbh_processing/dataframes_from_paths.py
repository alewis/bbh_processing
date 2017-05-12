import central_file
import parse_params
from os.path import expanduser, isfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
"""Construct a dataframe of initial eccentricity data.
"""

def read_old(path, df_to_concat):
    if path==None:
        return df_to_concat
    else:
        expanded_path = expanduser(path)
        if isfile(expanded_path):
            print "Appending data from", expanded_path
            old_dataframe = pd.read_pickle(expanded_path)
            new_dataframe = pd.concat([old_dataframe, df_to_concat], ignore_index=True)
            return new_dataframe
        else:
            return df_to_concat

def pickle_dataframe(
    func,
    datapaths = "~/.bbh/.datapaths", #file storing paths to the runs
    outpath = "~/.bbh/bbh_data.pkl", #where we will save the output
    oldpath1 = "~/.bbh/bbh_data_guil.pkl",
    oldpath2 = "~/.bbh/bbh_data_scinet.pkl"
    ):

    """ This runs the function func on all the data pointed to by datapaths.
        It then pickles the results at outpath.
    """
    #Read in the data paths
    expanded_data_path = expanduser(datapaths).strip()
    expanded_outpath = expanduser(outpath).strip()

    with open(expanded_data_path, 'r') as f:
      datapathnames = f.read().splitlines()
    for path in datapathnames:
      path = expanduser(path).strip()
    new_dataframe = func(datapathnames)
    new_dataframe = read_old(oldpath1, new_dataframe)
    new_dataframe = read_old(oldpath2, new_dataframe)
    new_dataframe.to_pickle(expanded_outpath)
    return new_dataframe

def checkfile(path):
    """Returns true if a file exists in path. Prints a warning otherwise.
       TODO: have this generate files if necessary and possible.
    """
    file_exists = isfile(path)
    if not file_exists:
        print "Warning - ", path, " not found."
    return file_exists

def orbital_dataframes(paths=None):
    """ Constructs, or loads, the dataframe.
    """
    lst = []
    for path in paths:
        omegapath = path+"Omega.pkl"
        parampath = path+"Params.input"
        hpath = path+"Horizons.h5"
        allfiles = True
        allfiles = checkfile(parampath)
        allfiles = checkfile(hpath)
        if allfiles:
            print "Adding data from ", path
            lst.append(parse_params.parse_params(parampath, 
                                                 hpath,
                                                 omegapath
                                                 )
                      )
    the_df = pd.DataFrame(lst)
    return the_df

if __name__ == "__main__":
    df = pickle_dataframe(orbital_dataframes)
