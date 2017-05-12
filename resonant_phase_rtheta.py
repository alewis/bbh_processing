#This creates plots of the 'resonant phase' of omega_theta vs. omega_r; defined
#as cos(int[0->t]{n*omega_theta + k*omega_r}dt). The argument to the cosine
#attains a maximum at the resonant time when the frequencies multipled by the
#resonant orders are equal (so the argument vanishes there). It also creates
#trajectory plots in r/theta space 1 radian in each direction around the resonant
#time.

import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import bbh_processing.resonant_phase as res
import bbh_processing.omega_theta_methods as omth
import bbh_processing.omega_r_methods as omr
import bbh_processing.hfile_tools as hf
import bbh_processing.angle_variables as av





if __name__=="__main__":
  #PRELIMINARY EXTRACTION OF DATA 
  pname = "Omega.pkl"
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)
  om_th = omth.omega_theta_from_fits(width=4)
  #om_r = omr.omega_r_from_ranges(OmegaMag,2)
  #om_r = omr.omega_r_from_ranges(OmegaMag,3)
  om_r = omr.omega_r_from_ranges(OmegaMag,4)
  #plt.show()
  #plt.plot(om_th[:,0],om_th[:,1]) 
  #plt.plot(om_r[:,0],om_r[:,1]) 
  #plt.show()
  hname = "Horizons.h5"
  time = hf.t(hname)
  
  theta = hf.theta_r_chi(hname)
  theta_arr = np.zeros((len(time),2))
  theta_arr[:,0] = np.copy(time)
  theta_arr[:,1] = np.copy(theta)
  
  rmag = hf.r_mag(hname)
  r_arr = np.zeros((len(time),2))
  r_arr[:,0] = np.copy(time)
  r_arr[:,1] = np.copy(rmag)
  
  #ratios = [(-2,3),(-1,4),(-1,5),(-3,4),(-2,5),(-1,6),(-4,5),(-5,6)]
  ratios = [(-2,3),(-3,4),(-4,5),(-5,6)]
  #ratios = [(-2,3),(-3,4),(-4,5),(-5,6)] 
  
  #COMPUTING THE RESONANT PHASES
  phases = [res.resonant_phase(om_th,om_r,r1,r2,0,len(om_r)) 
      for r1,r2 in ratios]  
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
  segs=[]
  for p,t_res in zip(cos_phases,tress):
    t1, tf = res.time_on_resonance_xrad(p,t_res,x=0.5) 
    segs.append(av.torus_plot_around_resonance_th_r(t1,tf))
    trajs.append(res.trajectory_around_resonance_xrad(theta_arr,r_arr,p,t_res,x=0.5))
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
  ax1w = plt.subplot2grid(grd,(0,3)) #angle
  ax1t = plt.subplot2grid(grd,(0,2)) #traj
  ax2w = plt.subplot2grid(grd,(1,3))
  ax2t = plt.subplot2grid(grd,(1,2))
  ax3w = plt.subplot2grid(grd,(2,3))
  ax3t = plt.subplot2grid(grd,(2,2))
  ax4w = plt.subplot2grid(grd,(3,3))
  ax4t = plt.subplot2grid(grd,(3,2))
  
    
  #format axes - phase plots
  for title,axp in zip(titles,axps):
    axp.set_ylabel(title+"\n\n"+r"$\cos(\phi_{res})$")
    #axp.set_title(title)
  ax4.set_xlabel(r"$t$") 
  axws = [ax1w,ax2w,ax3w,ax4w] 
  #format axes - angle plots
  for title,axw in zip(titles,axws):
    axw.set_xlim([0.0,2.0*np.pi]) 
    axw.set_ylim([0.0,2.0*np.pi]) 
    axw.set_xlabel(r"$W_r$")
    axw.set_ylabel(r"$W_{\theta}$"+"\n")
    axw.set_aspect('equal')
    #axw.set_title(title)
  
  axts = [ax1t,ax2t,ax3t,ax4t] 
  #format axes - trajectory plots
  for title,axt in zip(titles,axts):
    #axt.set_xlim([0.0,2.0*np.pi]) 
    #axt.set_ylim([0.0,2.0*np.pi]) 
    axt.set_xlabel(r"$r$")
    axt.set_ylabel(r"$\theta$")
    #axt.set_aspect('equal')
    #axt.set_title(title)
  
  bundle = zip(cos_phases,segs,trajs,tress,trngs,axps,axws,axts)
  ###MAKING THE PLOT
  for phase,segment,traj,tres,trng,axp,axw,axt in bundle:
    ti = trng[0]
    tf = trng[1]
    #the phase plots
    axp.plot(phase[:,0],phase[:,1]) 
    
    #lines around the resonance
    axp.plot((tres,tres),(-1,1),color="green")
    axp.plot((ti,ti),(-1,1),color="red")
    axp.plot((tf,tf),(-1,1),color="red")
    #grid of trajectory plots
    for s in segment:
      #the 'segment' is a bunch of lines segments, not one object, which is
      #why we loop here
      axw.plot(s[:,1],s[:,0],color="blue") 
    axt.plot(traj[:,2],traj[:,1],color="blue") 
  fig=plt.gcf()
  fig.set_size_inches(30,15)
  plt.savefig('image.pdf',format='pdf')  
  #plt.show()
