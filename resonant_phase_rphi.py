import numpy as np
from scipy import integrate
from scipy import interpolate as interp
import matplotlib.pyplot as plt
import bbh_processing.resonant_phase as res
import bbh_processing.omega_phi_methods as omphi
import bbh_processing.omega_r_methods as omr
import cPickle as pickle 
  

if __name__=="__main__":
  
  pname = "Omega.pkl"
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)

  #COMPUTING THE RESONANT PHASES
  phases = [res.resonant_phase(om_th,om_r,r1,r2,0,len(om_r)) for r1,r2 in ratios]  
  cos_phases = [] 
  for phase in phases:
    cphase = np.copy(phase)
    cphase[:,1] = np.cos(cphase[:,1])
    cos_phases.append(cphase) 
 
  #COMPUTING THE RESONANT TIMES AND TRAJECTORIES
  tress = [res.t_resonance(p) for p in phases] 
  titles = [str(r2)+":"+str(-r1) for r1,r2 in ratios]
  trngs = []
  trajs=[]
  for p,t_res in zip(cos_phases,tress):
    trajs.append(res.trajectory_around_resonance_1rad(phi_arr,r_arr,p,t_res))
    t1, tf = res.time_on_resonance_1rad(p,t_res) 
    trngs.append((t1,tf))
  
  
  ###SETTING UP THE PLOT
  from matplotlib import gridspec
  fig = plt.figure()
  #figtitle = "ecc0.3, inc 80"
  #make the grid
  grd = (4,4)
 
 
  #the phase plots
  ax1 = plt.subplot2grid(grd,(0,0),colspan=2)
  ax2 = plt.subplot2grid(grd,(1,0),colspan=2)#,sharex=ax1)
  ax3 = plt.subplot2grid(grd,(2,0),colspan=2)#,sharex=ax1)
  ax4 = plt.subplot2grid(grd,(3,0),colspan=2)#,sharex=ax1)
  axps = [ax1,ax2,ax3,ax4] 
  plt.setp(ax1,xticklabels=[]) 
  plt.setp(ax2,xticklabels=[]) 
  plt.setp(ax3,xticklabels=[]) 
  #the trajectory plots
  ax1t = plt.subplot2grid(grd,(0,2),rowspan=2)
  ax2t = plt.subplot2grid(grd,(0,3),rowspan=2)
  ax3t = plt.subplot2grid(grd,(2,2),rowspan=2)
  ax4t = plt.subplot2grid(grd,(2,3),rowspan=2)
  axts = [ax1t,ax2t,ax3t,ax4t] 
  
  bundle = zip(cos_phases,titles,trajs,tress,trngs,axps,axts)
  ###MAKING THE PLOT
  for phase,title,traj,tres,trng,axp,axt in bundle:
    ti = trng[0]
    tf = trng[1]
    axp.plot(phase[:,0],phase[:,1]) 
    axp.plot((tres,tres),(-1,1),color="green")
    axp.plot((ti,ti),(-1,1),color="red")
    axp.plot((tf,tf),(-1,1),color="red")
    axp.set_title(title)
    #grid of trajectory plots
    axt.plot(traj[:,2],traj[:,1]) 
    axt.set_ylabel("phi")
    axt.set_xlabel("r")
    axt.set_title(title)
  plt.show()
