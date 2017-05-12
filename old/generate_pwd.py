from matlabinterp import interp
import numpy as np
import itertools
import multiprocessing as mproc
import HorizonsTools  
import QuasiPeriodicTimeSeries as TS
import cPickle as pickle
import matplotlib.pyplot as plt
import scipy.signal
#slice = 5000

#t =np.asarray(HorizonsTools.Times())[4000:slice]
#r,rmag=HorizonsTools.SeparationVector()
#x = np.asarray(rmag)[4000:slice]
t = np.linspace(0,10*np.pi,5000)
x = scipy.signal.chirp(t,0,10,5000)
fftl = 2*len(x)+1
dt = t[1] - t[0]
freq = np.fft.fftfreq(fftl,dt)

if(t.shape[-1] != x.shape[-1]):
  raise ValueError(
      "Time and data should have same length.")
#x: input time series
#y: length of FFT to use

#get the length of the signal and the frequency vector
L = len(x)

#y is going to hold the (2D) output

#padding to extend the data to the fftl length
zp = fftl - 2*L + 1 
if zp < 0:
  raise ValueError(
      "fftl needs to be at least the length of the data")

#make sure we have data at the points we need
xinter = interp(x,2)
xinter = xinter[0:2*L-1] #why do we drop the last data point?



def make_one_pwd_column(k):
    mlk = k+1 #the index matlab uses
    print "k=", k, " L-k =", L-k
    
    #Make the x(t+tau/2) for the Wigner integral.
    if k<L/2:
      xx1 = np.zeros(2*L-1)
      #Take a slice from the 0th to the 2(L/2+mlk)-1th entry of x. 
      padlength = 2*(L/2+mlk)-1 
      datalength = 2*(L/2+mlk)-1
      xx1[0:datalength] = xinter[0:datalength]
      
    else:
      xx1 = np.zeros(2*L-1) 
      data_start = 2*(-L/2+mlk)-2
      data_end = 2*L-1
      xx1[0:(data_end-data_start)] = xinter[data_start:data_end]
    
    #Lin alg. trickery to get the x(t-tau/2)*
    xx2 = np.conj(xx1[::-1]) #[::-1] creates a view to the reversed array
    
    #Multiply them together elementwise: x(t+tau/2) * x(t-tau/2)*
    xx = xx1*xx2 
    Lxx = xx[L-1:2*L-1] 
    xxR = xx[0:L-1] 
    N = len(Lxx)+zp+len(xxR)
    xxzp = np.zeros(N)
    xxzp[:len(Lxx)]=Lxx
    xxzp[len(Lxx)+zp:]=xxR
    return (k,np.abs(np.fft.fft(xxzp).transpose()))
    #return np.fft.fft(xxzp).transpose()

def start_proc():
  print "Starting: ", mproc.current_process().name

#Handle the parallelization
def main():
  pool_size = mproc.cpu_count()*2 
  pool_size = 1
  print pool_size
  pool = mproc.Pool(processes = pool_size,
                    initializer = start_proc,
                    maxtasksperchild=2,
                    )
  procids = list(range(L))
  pool_outs=pool.map(make_one_pwd_column,procids)
  pool.close()
  pool.join()
  
  y = np.zeros((fftl,L))
  for slice in pool_outs:
    y[:,slice[0]] = slice[1]
  
  plt.imshow(y) 
  plt.show()
  #g = open("PWD.pkl",'w')
  #pickle.dump([t,freq,y],g)
  #g.close()

if __name__ == '__main__':
  #generate data
  main()
  print "Done!"
    
