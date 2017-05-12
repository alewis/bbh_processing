#Various functions useful for analyizing black hole data
import h5py
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
#######################################
#Deal with W files
#######################################
def omega_phi_from_w(wname="W1.6.dat"):
  return np.loadtxt(wname,usecols=(0,1))

def omega_r_from_w(wname="W1.6.dat"):
  return np.loadtxt(wname,usecols=(0,2))

def k_from_w(wname="W1.6.dat"):
  return np.loadtxt(wname,usecols=(0,3))





#def interp_and_apply_function(small_time,big_time,small,big,function,
#    order="small first"):
#  big_made_small = interp_onto(big,big_time,small_time,'cubic')
#  if "small first"==order:
#    return function(small,big_made_small)
#  else if "big first"==order:
#    return function(big_made_small,small)
#  else:
#    raise ValueError("keyword 'order' must be 'small first' or 'big first'")

#STUFF FROM FLUX FILES

def FluxTime(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(0,))
def dEdt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(1,))
def dJxdt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(2,))
def dJydt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(3,))
def dJzdt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(4,))
def dPxdt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(5,))
def dPydt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(6,))
def dPzdt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(7,))
def E(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(8,))
def Jx(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(9,))
def Jy(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(10,))
def Jz(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(11,))
def Px(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(12,))
def Py(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(13,))
def Pz(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(14,))
  
def Jdot(fname="Fluxes.dat"):
  return np.vstack((dJxdt(fname),dJydt(fname),dJzdt(fname))).transpose()
def J(fname="Fluxes.dat"):
  return np.vstack((Jx(fname),Jy(fname),Jz(fname))).transpose()
def Pdot(fname="Fluxes.dat"):
  return np.vstack((dJxdt(fname),dJydt(fname),dJzdt(fname))).transpose()
def P(fname="Fluxes.dat"):
  return np.vstack((Px(fname),Py(fname),Pz(fname))).transpose()



  
def sinewavegueses(data):
  turns,turn_vals,rngs = turning_points(data)
  #amplitude guesses - half the peak to trough distance
  zipvals = zip(turn_vals[0::1],turn_vals[1::1])
  amps = [math.fabs((z[0]-z[1])/2) for z in zipvals]
  
  #offset guesses - the mean of the data in each range
  offsets = [np.mean(data[r[0]:r[1]]) for r in rngs]
  
  #frequency guesses - 1/ twice the peak to trough time
  freqs = [1./math.fabs(2*(z[0]-z[1])) for z in zip(turns[0::1],turns[1::1])]
  
  #phase guesses - pi/2 for a peak, 3pi/2 for a trough
  phases=[]
  pphase = np.pi/2
  tphase = 3*np.pi/2
  end = len(zipvals)
  for i in range(0,end-1):
    if zipvals[i] > zipvals[i+1]:
      phases.append(pphase*(1/freqs[i]))
    else:
      phases.append(tphase*(1/freqs[i]))
  if zipvals[end-1] < zipvals[end-2]:
      phases.append(pphase*(1/freqs[end-1]))
  else:
      phases.append(tphase*(1/freqs[end-1]))
  
  return amps,offsets,freqs,phases 

def AxialCMMotion(q,fname="Horizons.h5",mtot=1):
  f=h5py.File(fname,'r')
  m1 = mtot*(1-1/(q+1))
  m2 = mtot/(q+q)
  xA = np.array(f["/AhA.dir"]["CoordCenterInertial.dat"])
  xB = np.array(f["/AhB.dir"]["CoordCenterInertial.dat"])
  #coordinate components
  Xa = m1*xA[:,1]/mtot
  Ya = m1*xA[:,2]/mtot
  Za = m1*xA[:,3]/mtot
  Xb = m2*xB[:,1]/mtot
  Yb = m2*xB[:,2]/mtot
  Zb = m2*xB[:,3]/mtot
  
  xcm = Xa-Xb
  ycm = Ya-Yb
  zcm = Za-Zb
  #Xa -= xcm[0]
  #Xa -= ycm[0]
  #Xa -= zcm[0]
  plt.plot(xA[:,0],zcm, label="z")
  
  plt.plot(xA[:,0],np.sqrt(xcm**2+ycm**2),label="rho")
  plt.legend()
  plt.show()

def MonotonicFit(data,time):
  import scipy.optimize as opt
  turns,turn_vals,rngs=turning_points(data)
  amps,offsets,freqs,phases = sinewavegueses(data)
  #plt.plot(time,data)
  print len(turns), len(offsets)
  rngtime=[]
  for r in rngs:
    rngtime.append(time[(r[0]+r[1])/2])
  
  powerlaw = lambda p, x: p[0]*(p[1]-x)**p[2]+p[3]
  powererror = lambda p, x, y: powerlaw(p,x) - y
  
  #p0=[1,1,0.25]
  #print p1
  #plt.plot(rngtime,offsets)
  pwrp0 = [-100,time[-1],-0.25,offsets[0]]
  
  pwrp1, success = opt.leastsq(powererror,pwrp0[:],args=(rngtime,offsets))
  smoothoffsets=[powerlaw(pwrp1,t) for t in time]
  #plt.plot(time,smoothoffsets)
  #plt.show()
  #now interpolate
  
  return np.vstack((time,smoothoffsets))

#rolling sinusoidal fit over time,data timeseries
def RollingPowerLawSinusoid(data,time):
  import scipy.optimize as opt
  turns,turn_vals,rngs=turning_points(data)
  amps,offsets,freqs,phases = sinewavegueses(data)
  plt.plot(time,data)
  print len(turns), len(offsets)
  rngtime=[]
  for r in rngs:
    r*rtnerAnbdngtime.append(time[(r[0]+r[1])/2])
  
  powerlaw = lambda p, x: p[0]*(p[1]-x)**p[2]+p[3]
  powererror = lambda p, x, y: powerlaw(p,x) - y
  
  #p0=[1,1,0.25]
  #print p1
  #plt.plot(rngtime,offsets)
  pwrp0 = [-100,time[-1],-0.25,offsets[0]]
  
  pwrp1, success = opt.leastsq(powererror,pwrp0[:],args=(rngtime,offsets))
  offsmth = lambda x: pwrp1[0]*(pwrp1[1]-x)**pwrp1[2] + pwrp1[3] 
  #smoothoffsets=[powerlaw(p1[:],t) for t in rngtime]
  #plt.show()
  #now fit sinusoids over about two oscillation periods
  #need to fit amplitude AMP, offset OFF, and phase PHA
  #to do this we need 'reasonable first guesses' of each

  #first guess of OFF is the data set's mean; this should
  #be stable enough that we only need to do this once
  pi = math.pi
  fittimes=[]
  fit=[]
  om_t = []
  om_t_times=[]
  diff = []
  
  fitfunc = lambda p, x: offsmth(x)+(p[0]*np.sin(2*pi*p[1]*x + p[2])) 
  errfunc = lambda p, x, y: fitfunc(p,x) - y
  #we're at a peak - do the initial fit
  start = rngs[0][0]
  end = rngs[len(rngs)-1][1]
  b = rngs[0][1] 
  A0 = amps[0]
  B0 = freqs[0]
  C0 = phases[0]
  #D0 = powerlaw(pwrp1[:],time[start])
  p0 = [A0, B0, C0]
  for j in range(start,end):
    #print range
    om_t_times.append(time[j])
    p1, success = opt.leastsq(errfunc, p0[:], args=(time[j:b],data[j:b]))
    print "t: ", time[j], "A1: ", p1[0], "    B1: " , p1[1], "    C1: ", p1[2]
    p0=p1  
    
    fit.append(fitfunc(p1[:],time[j]))
    om_t.append(p1[1])
    b+=1
  return om_t_times,om_t,fit    

#rolling sinusoidal fit over time,data timeseries
def RollingSinusoid(data,time):
  import scipy.optimize as opt
  turns,turn_vals,rngs=turning_points(data)
  amps,offsets,freqs,phases = sinewavegueses(data)
  #RollingPowerLawSinusoid(data,time)
#now fit sinusoids over about two oscillation periods
#need to fit amplitude AMP, offset OFF, and phase PHA
#to do this we need 'reasonable first guesses' of each

#first guess of OFF is the data set's mean; this should
#be stable enough that we only need to do this once
  pi = math.pi
  fittimes=[]
  fit=[]
  om_t = []
  om_t_times=[]
  diff = []
  
  fitfunc = lambda p, x: p[0]*np.sin(2*pi*p[1]*x + p[2])+np.mean(offsets)
  errfunc = lambda p, x, y: fitfunc(p,x) - y
  #we're at a peak - do the initial fit
  start = rngs[0][0]
  end = rngs[len(rngs)-1][1]
  b = rngs[0][1] 
  A0 = amps[0]
  B0 = freqs[0]
  C0 = phases[0]
  p0 = [A0, B0, C0]
  for j in range(start,end):
    #print range
    om_t_times.append(time[j])
    p1, success = opt.leastsq(errfunc, p0[:], args=(time[j:b],data[j:b]))
    print "t: ", time[j], "A1: ", p1[0], "    B1: " , p1[1], "    C1: ", p1[2]
    p0=p1  
    
    fit.append(fitfunc(p1[:],time[j]))
    om_t.append(p1[1])
    b+=1
  return om_t_times,om_t,fit    

#rolling sinusoidal fit over time,data timeseries
def RollingSinusoid(data,time):
  import scipy.optimize as opt
  turns,turn_vals,rngs=turning_points(data)
  amps,offsets,freqs,phases = sinewavegueses(data)
  #RollingPowerLawSinusoid(data,time)
#now fit sinusoids over about two oscillation periods
#need to fit amplitude AMP, offset OFF, and phase PHA
#to do this we need 'reasonable first guesses' of each

#first guess of OFF is the data set's mean; this should
#be stable enough that we only need to do this once
  pi = math.pi
  fittimes=[]
  fit=[]
  om_t = []
  om_t_times=[]
  diff = []
  
  fitfunc = lambda p, x: p[0]*np.sin(2*pi*p[1]*x + p[2])+np.mean(offsets)
  errfunc = lambda p, x, y: fitfunc(p,x) - y
  #we're at a peak - do the initial fit
  start = rngs[0][0]
  end = rngs[len(rngs)-1][1]
  b = rngs[0][1] 
  A0 = amps[0]
  B0 = freqs[0]
  C0 = phases[0]
  p0 = [A0, B0, C0]
  for j in range(start,end):
    #print range
    om_t_times.append(time[j])
    p1, success = opt.leastsq(errfunc, p0[:], args=(time[j:b],data[j:b]))
    print "t: ", time[j], "A1: ", p1[0], "    B1: " , p1[1], "    C1: ", p1[2]
    p0=p1  
    
    fit.append(fitfunc(p1[:],time[j]))
    om_t.append(p1[1])
    b+=1
  return om_t_times,om_t,fit    

