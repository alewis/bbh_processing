import numpy as np
def make_one_pwd_column(k, xinter, L, zp,y):
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
    
    y[:,k]= np.fft.fft(xxzp).transpose()
