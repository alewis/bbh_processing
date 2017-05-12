from make_one_pwd_column import make_one_pwd_column

def PseudoWigner(t,x, fftl):
  from matlabinterp import interp
  if(t.shape[-1] != x.shape[-1]):
    raise ValueError(
        "Time and data should have same length.")
  #x: input time series
  #y: length of FFT to use

  #get the length of the signal and the frequency vector
  L = len(x)

  #y is going to hold the (2D) output
  y = np.zeros((fftl,L))
  #xx is going to hold the shifted input for the Wigner filter
  xxlen = 2*L+1
  xx = np.zeros(xxlen)
  
  #padding to extend the data to the fftl length
  zp = fftl - 2*L + 1 
  if zp < 0:
    raise ValueError(
        "fftl needs to be at least the length of the data")

  #make sure we have data at the points we need
  xinter = interp(x,2)
  xinter = xinter[0:2*L-1] #why do we drop the last data point?
  
  for i in range(L):
    y[:,i] = make_one_pwd_column(i)
  dt = t[1] - t[0]
  freq = np.fft.fftfreq(fftl,dt)
  return freq,y
