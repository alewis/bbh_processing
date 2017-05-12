import cPickle as pickle
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from peak_finders import turning_points, monotonic_peaks,monotonic_peaks2
def ecc(ra,rp):
  return np.fabs((ra-rp))/(ra+rp)

def load_ecc_data(sepname,pname,mode="coord"):
  #Load in orbital frequency data 
  f = open(pname)
  Omega,OmegaMag=pickle.load(f)
  f.close()
  cols=()

  #Load in proper or coordinate separation depending on choice of mode 
  if mode=="proper":
    cols = (0,1)
  elif mode=="coord":
    cols = (0,2)
  else:
    raise ValueError("Mode must be 'coord' or 'proper'")
  t, r = np.loadtxt(sepname,usecols=cols,unpack=True)
  return t,r,OmegaMag

#Generate a list of eccentricities at the apastrons of a data set,
#and store them along with the associated orbital frequencies in a file
def make_ecc_list(
                   sepname,
                   pname,
                   mode    = "coord"
                  ):
 
  t,r,OmegaMag = load_ecc_data(sepname,pname,mode)
  
  #trim the noisy bit off the beginning
  t = t[200:]
  r = r[200:]
  turns, turn_vals, rngs = turning_points(r,first_left=0,first_right=1,stride=1)

  #trim so we start at apastron
  try:
    if (turn_vals[0]<turn_vals[1]):
      turns = turns[1:]
      turn_vals = turn_vals[1:]
      rngs = rngs[1:]
  except IndexError:
    print "Couldn't find enough peaks; breaking..."
    return dict()

  #and finish at periastron
  if (turn_vals[-2]<turn_vals[-1]):
    turns = turns[:-2]
    turn_vals=turn_vals[-2]
    rngs=rngs[-2]
  ecclist = []
  rlist = []
  freqlist = []
  for i in range(0,len(turns),2):
    try: 
      ecclist.append(ecc(r[turns[i]],r[turns[i+1]]))
      freqlist.append(OmegaMag[:,1][turns[i]])
      rlist.append(r[turns[i]]) 
    except IndexError:
      print "Warning: index was out of bounds."
  return rlist,freqlist,ecclist


#Returns a numpy array mapping (r,freq) pairs to ecc
def make_ecc_array(
                   sepname,
                   pname,
                   coordmode    = "coord",
                   eccmode      = "interp"
                  ):
  if (eccmode=="old"):
    rlist,freqlist,ecclist = make_ecc_list(sepname,pname,coordmode)
  elif (eccmode=="interp"):
    errcode,tshort,rlist,freqlist,ecclist = ecc_from_interpolated_peaks2(
        sepname,pname,coordmode) 
    if (errcode==1): return errcode,0 
  elif (eccmode=="peak"):
    errcode,tshort,rlist,freqlist,ecclist = ecc_peaks(sepname,pname,coordmode) 
    if (errcode==1): return errcode,0 
  elif (eccmode=="interp_app"):
    errcode,tshort,rlist,freqlist,ecclist = ecc_peaks(sepname,pname,coordmode,True) 
    if (errcode==1): return errcode,0 
  the_array = np.asarray(zip(tshort,rlist,freqlist,ecclist))
  print the_array
  return 0, the_array


#Returns a dictionary mapping (r,freq) pairs to ecc
def make_ecc_dict(
                   sepname,
                   pname,
                   coordmode    = "coord",
                   eccmode      = "interp"
                  ):
  if (eccmode=="old"):
    rlist,freqlist,ecclist = make_ecc_list(sepname,pname,coordmode)
  elif (eccmode=="interp"):
    try:
      tshort,rlist,freqlist,ecclist = ecc_from_interpolated_peaks2(sepname,pname,coordmode) 
    except:
      print "INTERPOLATION FAILED"
      return dict()
  #zip r and freq into duples to use as keys 
  keys = zip(rlist,freqlist)
  #and return the dictionary
  return dict(zip(keys,ecclist))

#Return the eccentricity at the peaks of the time series along with the times,
#separations, and orbital frequencies at those peaks. The eccentricity is estimated
#by doing a spline interpolation between the peaks and the troughs of the data 
#set and then computing (ra-rp)/(ra+rp) from the results.
def ecc_peaks(sepname,pname,mode="coord",app_only=False):
  #read the data from the separation (HorizonSepMeasures.dat) and orbital frequency
  #(Omega.pkl) files. 'Mode' chooses between coordinate and proper separation
  #('coord' or 'proper').
  t,r,OmegaMag = load_ecc_data(sepname,pname,mode)
 
  #errcode is 1 if there weren't enough peaks to make the splines.
  #left and right bound the range where both splines are simultaneously
  #interpolated.
  #peakspline interpolates the peaks; trougspline the troughs
  errcode,left,right,peakspline,troughspline = monotonic_peaks2(t,r)
  if errcode==1:
    return errcode,0,0,0,0
  
  
  #The indices and values of the peaks.
  turns,turn_vals = turning_points(r) 
  if(app_only==True): 
    vals1 = sum(turn_vals[::2])
    vals2 = sum(turn_vals[1::2])
    #Apastron will have greater summed separations 
    if(vals1>vals2):
      turns=turns[::2]
      turn_vals=turn_vals[::2]
    else:
      turns=turns[1::2]
      turn_vals=turn_vals[1::2]
  #Trim them so that they have the same range as OmegaMag
  #and the spline ranges
  tarr = np.asarray(turns)
  varr = np.asarray(turn_vals)
  
  om_idx = len(OmegaMag)
  #choose the rightmost good index
  rcut = min(len(OmegaMag),(right-1))
  #prepare an array of 2 length = len(tarr) columns of bools 
  trimbool = np.ones((len(tarr),2),dtype=bool)
  #column 0 is true for indices >= left 
  trimbool[:,0] = tarr>= left
  #column 1 is true for indices <= rcut
  trimbool[:,1] = tarr<= rcut
  
  #np.all does an array-wise 'and'
  trimall = np.all(trimbool,axis=1)
  
  #now we actually trim
  
  #trimall is the list of bools determining whether a peak is within the
  #range or not, and is the appropriate filter for quantities stored at peaks
  #(turns, turn_vals, and derivatives). ttrim is the indices of the peaks
  #after trimming, and is appropriate for full time series
  ttrim = tarr[trimall]
  peakt = t[ttrim].tolist() 
  peakr = varr[trimall].tolist()
  peakom = OmegaMag[ttrim,1].tolist()
  peakecc = [ecc(peakspline(ti), troughspline(ti)) for ti in peakt]
  #plt.plot(peakt,peakr)
  #plt.plot(t,r)
  #plt.show()
  return errcode, peakt,peakr,peakom,peakecc

#As eccpeaks, but don't return the frequencies and feed in the separation data
#directly (instead of a file name)
def ecc_peaks_from_data(t,r):
  #read the data from the separation (HorizonSepMeasures.dat) and orbital frequency
  #(Omega.pkl) files. 'Mode' chooses between coordinate and proper separation
  #('coord' or 'proper').
 
  #errcode is 1 if there weren't enough peaks to make the splines.
  #left and right bound the range where both splines are simultaneously
  #interpolated.
  #peakspline interpolates the peaks; trougspline the troughs
  errcode,left,right,peakspline,troughspline = monotonic_peaks2(t,r)
  if errcode==1:
    return errcode,0,0,0,0
  
  
  #The indices and values of the peaks.
  turns,turn_vals = turning_points(r) 
  
  #select only the apastron peaks 
  vals1 = sum(turn_vals[::2])
  vals2 = sum(turn_vals[1::2])
  #Apastron will have greater summed separations 
  if(vals1>vals2):
    turns=turns[::2]
    turn_vals=turn_vals[::2]
  else:
    turns=turns[1::2]
    turn_vals=turn_vals[1::2]
  #Trim them so that they have the same range as the spline ranges
  tarr = np.asarray(turns)
  
  #choose the rightmost good index
  rcut = right-1
  #prepare an array of 2 length = len(tarr) columns of bools 
  trimbool = np.ones((len(tarr),2),dtype=bool)
  #column 0 is true for indices >= left 
  trimbool[:,0] = tarr>= left
  #column 1 is true for indices <= rcut
  trimbool[:,1] = tarr<= rcut
  
  #np.all does an array-wise 'and'
  trimall = np.all(trimbool,axis=1)
  
  #now we actually trim
  
  #trimall is the list of bools determining whether a peak is within the
  #range or not, and is the appropriate filter for quantities stored at peaks
  #(turns, turn_vals, and derivatives). ttrim is the indices of the peaks
  #after trimming, and is appropriate for full time series
  ttrim = tarr[trimall]
  peakt = t[ttrim].tolist() 
  peakecc = [ecc(peakspline(ti), troughspline(ti)) for ti in peakt]
  return errcode, peakt,peakecc

def ecc_from_interpolated_peaksdata(t,r):
  errcode,left,right,peakspline,troughspline = monotonic_peaks2(t,r)
  if errcode==1:
    return errcode,0,0,0,0
  else: 
    tshort = t[left:right]
    ecc_list = ecc(peakspline(tshort),troughspline(tshort))
    plt.plot(tshort,peakspline(tshort))
    plt.plot(tshort,troughspline(tshort))
    print ecc_list
    return errcode,tshort,ecc_list

def ecc_from_interpolated_peaks2(sepname,pname,mode="coord"):
  t,r,OmegaMag = load_ecc_data(sepname,pname,mode)
  errcode,left,right,peakspline,troughspline = monotonic_peaks2(t,r)
  if errcode==1:
    return errcode,0,0,0,0
  else: 
    tshort = t[left:right]
    ecc_list = ecc(peakspline(tshort),troughspline(tshort))
    print ecc_list
    rshort = r[left:right]
    omshort = OmegaMag[left:right,1]
    return errcode,tshort,rshort,omshort,ecc_list


def ecc_from_interpolated_peaks(sepname,pname,mode="coord"):
  t,r,OmegaMag = load_ecc_data(sepname,pname,mode)
  leftrange,rightrange,peakspline,troughspline = monotonic_peaks(t,r)
  #convert leftrange,rightrange to indices
  idxleft = (np.abs(t-leftrange)).argmin()
  idxright = (np.abs(t-rightrange)).argmin()

  tshort = t[idxleft:idxright]
  ecc_list = ecc(peakspline(tshort),troughspline(tshort))
  rshort = r[idxleft:idxright]
  omshort = OmegaMag[idxleft:idxright,1]
  #short_t = t[leftrange:rightrange]
  #ecc_list = ecc(peakspline(short_t),troughspline(short_t))
  #r_int = interp_onto(r,t,short_t)
  #om_int = interp_onto(OmegaMag[:,1],OmegaMag[:,0],short_t)

  #rlist = [return_interpolated_value(t,st,50,r) for st in short_t]
  #omlist =[return_interpolated_value(t,st,50,OmegaMag[:,1]) for st in short_t]
  return tshort,rshort,omshort,ecc_list
#tshort,rlist, omlist, ecclist = ecc_from_interpolated_peaks("HorizonSepMeasures.dat","Omega.pkl")
#plt.plot(tshort,rlist)
#plt.plot(tshort,ecclist)
#plt.show()


