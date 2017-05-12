import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as opt
import HorizonsTools  
import cPickle as pickle

theta=HorizonsTools.AngleWithSpin()
time =HorizonsTools.Times()
f=open("Omega.pkl",'r')
Omega,OmegaMag = pickle.load(f)
f.close()

#get all the turning points of the curve
turns = []
turn_vals = []
for pt in range(1,len(theta)-1):
    if (theta[pt] > theta[pt-1]):
        if (theta[pt] > theta[pt+1]):
            turns.append(pt)
            turn_vals.append(theta[pt])
    elif (theta[pt] < theta[pt-1]):
        if (theta[pt] < theta[pt+1]):
            turns.append(pt)
            turn_vals.append(theta[pt])

#get ranges of two periods each
rngs = zip(turns[1::4],turns[5::4])

#now fit sinusoids over about two oscillation periods
#need to fit amplitude AMP, offset OFF, and phase PHA
#to do this we need 'reasonable first guesses' of each

#first guess of OFF is the data set's mean; this should
#be stable enough that we only need to do this once
pi = math.pi

om_t = []
om_t_times=[]
diff = []
for range in rngs:
    a = range[0]
    b = range[1]
    #print range
    mid = a+(b-a)/2
    om_t_times.append(time[mid])
    fitfunc = lambda p, x: p[0]*np.sin(2*pi*p[1]*x + p[2]) + p[3]
    errfunc = lambda p, x, y: fitfunc(p,x) - y
    #y = A * sin(t*B + C) + D
    A0 =2*np.std(theta[a:b])/(2**0.5)
    period = (time[b]-time[a]) / 2.0
    
    B0 = 1.0/period
    
    C0 = 0
    D0 = np.mean(theta[a:b]) 
    #print "A0: ", A0, "    B0: " , B0, "     C0: ", C0, "    D0: " , D0

    p0 = [A0, B0, C0, D0]
    p1, success = opt.leastsq(errfunc, p0[:], args=(time[a:b],theta[a:b]))
    #print "A1: ", p1[0], "    B1: " , p1[1], "    C1: ", p1[2], "    D1: " , p1[3]
    om_t.append(p1[1])
    diff.append(OmegaMag[mid,1]/(2.0*pi) - p1[1])

#print OmegaMag[2,1]
 
#plt.ylabel("Ordinary frequency")
#plt.xlabel("Time (CM frame)")
#plt.title("Frequency of spin vector oscillations")
#plt.plot(OmegaMag[:,0],OmegaMag[:,1]/(2.0*pi), label = 'Orbital frequency')
#plt.plot(time,theta, label = 'theta')
#plt.plot(om_t_times,om_t, label = 'Spin frequency (Om_t)')
#plt.plot(om_t_times,diff, "r+", label = 'Difference')
#plt.legend()
#plt.show()
