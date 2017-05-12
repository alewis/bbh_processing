#omega_r: the frequency of the peak to peak oscillations in omega total
import numpy as np
import bbh_processing.utilities as bbh
import frequency_tools as freq


#Define omega_r as 1/ the time between successive peaks
def omega_r_from_ranges(OmegaMag, width=2):
  return freq.freq_from_ranges(OmegaMag,width=width) 

def omega_r_from_averaged_ranges(OmegaMag,minwidth=2,maxwidth=8):
  return bbh.mean_windowed_series(freq.freq_from_ranges,OmegaMag,
      minwidth,maxwidth,2)

def omega_r_from_r(t,r,minwidth=2,maxwidth=8):
  r_arr = np.zeros((len(t),2))
  r_arr[:,0] = t
  r_arr[:,1] = r
  return bbh.mean_windowed_series(freq.freq_from_ranges,r_arr,
      minwidth,maxwidth,2)
