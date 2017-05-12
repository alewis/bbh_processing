#return a plot of J projected onto S
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
import equatorial_freq_methods as eqfreq

def plot_krphvt(outdict, names, labels, colours, ax):
    for name, label, colour in zip(names, labels, colours):
        t_pk = outdict[name]["t_pk"]
        t_tr = outdict[name]["t_tr"]
        krph_pk = outdict[name]["krph_pk"]
        krph_tr = outdict[name]["krph_tr"]
        ax.plot(t_tr, krph_tr, marker="*", label=label+", tr", c=colour)
        ax.plot(t_pk, krph_pk, marker="*", label=label+", pk", c=colour,
                linestyle="--")
def plot_krphvph(outdict, names, labels, colours, ax):
    for name, label, colour in zip(names, labels, colours):
        omph_pk = outdict[name]["omph_pk"]
        omph_tr = outdict[name]["omph_tr"]
        krph_pk = outdict[name]["krph_pk"]
        krph_tr = outdict[name]["krph_tr"]
        ax.plot(omph_tr, krph_tr, marker="*", label=label+", tr", c=colour)
        ax.plot(omph_pk, krph_pk, marker="*", label=label+", pk", c=colour,
                linestyle="--")
        # ax.fill_between(t_pk, krph_pk, krph_tr)


#orbital frequency vector (L) dotted with chi
def L_dot_chi(hname="Horizons.h5", pname="Omega.pkl"):
  import cPickle as pickle
  import bbh_processing.hfile_tools as hf
  import bbh_processing.utilities as bbh
  f = open(pname)
  L,OmMag = pickle.load(f)
  f.close()
 
  time = hf.t(hname)
  chi = hf.chi_inertial(hname)
  #print time[10:-11] compensates for points lost to the omega finder
  #print L[:,0]
  ldotchi = bbh.timedot(L[:,1:],chi[10:-11,:])
  the_plot = plt.plot(L[:,0],ldotchi)
  plt.xlabel("Time(M)")
  plt.ylabel("L . chi")
  plt.title("L dot Chi")
  return the_plot

#orbital eccentricity via r_a r_p definition
def ecc_from_app_peri(hname="Horizons.h5"):
  import bbh_processing.hfile_tools as hf
  time, ecc = hf.ecc_from_app_peri(hname)
  the_plot = plt.plot(time,ecc)
  plt.xlabel("Time(M)")
  plt.ylabel("Eccentricity")
  plt.title("Eccentricity from r_a r_p")
  return the_plot

#orbital frequency
def omega_mag(pname="Omega.pkl"):
  import cPickle as pickle
  f = open(pname)
  Om,OmMag = pickle.load(f)
  f.close()
  the_plot = plt.plot(OmMag[:,0],OmMag[:,1])
  plt.xlabel("Time(M)")
  plt.ylabel("Orbital Frequency (1/M)")
  plt.title("Orbital Frequency")
  return the_plot

#theta frequency
def omega_theta_plot(pname="Omega.pkl", hname="Horizons.h5"):
  import cPickle as pickle
  from bbh_processing.combined_tools import omega_theta
  f = open(pname)
  om_th_t, om_th = omega_theta(pname,hname)
  f.close()
  the_plot = plt.plot(om_th_t,om_th)
  plt.xlabel("Time(M)")
  plt.ylabel("Theta frequency (1/M)")
  plt.title("Omega Theta")
  return the_plot

def omegas_from_abdul(wname="W1.6.dat", hname="Horizons.h5",pname="Omega.pkl"):
  import cPickle as pickle
  from bbh_processing.combined_tools import omega_theta
  import bbh_processing.wfile_tools as wf
  f = open(pname)
  Om,OmMag = pickle.load(f)
  f.close()
  wt = wf.t(wname)
  om_phi=wf.omega_phi(wname)
  om_r=wf.omega_r(wname)
  om_th_t, om_th = omega_theta(hname)
  
  
  #om_r_t, om_r = omega_r(pname)
  #om_phi_t, om_phi = omega_phi(pname)
  lab=['Om_th','Om_phi','Om_r']
  the_plot = plt.plot(om_th_t,om_th,wt,om_phi,wt,om_r,
      OmMag[:,0],OmMag[:,1])
  
  plt.legend(lab, loc='upper left',fontsize=8)
  plt.xlabel("Time(M)")
  plt.ylabel("Frequency (1/M)")
  plt.title("Omega")
  return the_plot


def omegas(pname="Omega.pkl", hname="Horizons.h5", equatorial=False):
  import cPickle as pickle
  import bbh_processing.combined_tools as ct 
  f = open(pname)
  Om,OmMag = pickle.load(f)
  f.close()
  
  lab=['Om_phi','Om_r']
  if equatorial==False:
    om_th_t, om_th = ct.omega_theta(hname)
    lab=['Om_th','Om_phi','Om_r']
  om_r_t, om_r = ct.omega_r(pname)
  #om_phi_t, om_phi = omega_phi(pname)
  om_phi_t,om_phi = ct.omega_phi(pname)
  
  if equatorial:
    the_plot = plt.plot(om_phi_t,om_phi,om_r_t,om_r,
      OmMag[:,0],OmMag[:,1])
  else:
    the_plot = plt.plot(om_th_t,om_th,om_phi_t,om_phi,om_r_t,om_r,
      OmMag[:,0],OmMag[:,1])
  
  plt.legend(lab, loc='upper left',fontsize=8)
  plt.xlabel("Time(M)")
  plt.ylabel("Frequency (1/M)")
  plt.title("Omega")
  return the_plot

def gen_ks(pname="Omega.pkl",hname="Horizons.h5",equatorial=False):
  import cPickle as pickle
  import bbh_processing.utilities as bbh
  import numpy as np
  from bbh_processing.combined_tools import omega_theta,omega_r,omega_phi
  f = open(pname)
  Om,OmMag = pickle.load(f)
  f.close()
  
  om_r_t, om_r = omega_r(pname)
  om_phi_t, om_phi = omega_phi(pname)
  k_r_phi = np.divide(om_r,om_phi)
  
  if equatorial:
    return om_phi, om_r_t, k_r_phi 

  else:
    #om_th has different shape from the others
    #(since it both starts and ends later) 
    #we therefore need to trim om_phi and om_r from the left to match the 
    #initial time of om_th, and om_th from the right to match the final time
    #of its counterparts
    
    om_th_t, om_th = omega_theta(hname)
    
    #first we get the values of the time extents we want
    initial_time = om_th_t[0]
    final_time = om_r_t[-1]
 
    #now the index of om_r and om_phi corresponding to the initial time
    initial_index_r = (np.where(om_r_t==initial_time)[0][0],-1)[0]
    initial_index_phi = (np.where(om_phi_t==initial_time)[0][0],-1)[0]
    #now the index of om_th corresponding to the final time
    final_index_th = (np.where(om_th_t==final_time)[0][0],-1)[0]
    om_r_trimmed = om_r[initial_index_r:-1]
    om_phi_trimmed = om_phi[initial_index_phi:-1]
    om_th_trimmed = om_th[0:final_index_th]

    k_r_th = np.divide(om_r_trimmed,om_th_trimmed)
    k_th_phi = np.divide(om_th_trimmed,om_phi_trimmed)
    trim_time = om_r_t[initial_index_r:-1]
    return (om_phi, om_phi_trimmed,om_r_t, k_r_phi, 
        trim_time, k_r_th, trim_time, k_th_phi) 


def get_resonant_times(pname="Omega.pkl",hname="Horizons.h5",
    equatorial=False,freq=False):
  import bbh_processing.utilities as bbh
  ratios = [1., 0.5, 2./3., 3./4., 4./5.] 
  r_phi_crossings = []
  if equatorial:
    om_phi, om_r_t, k_r_phi = gen_ks(pname,hname,equatorial) 
    for r in ratios:
      if freq:
        for x in bbh.invert(k_r_phi,om_phi,r):
          r_phi_crossings.append(x)
      else:
        for x in bbh.invert(k_r_phi,om_r_t,r):
          r_phi_crossings.append(x)
    return r_phi_crossings 
  else:
    r_th_crossings=[]
    th_phi_crossings=[]
    
    (om_phi,om_phi_trimmed,om_r_t, k_r_phi, trim_time, k_r_th, trim_time, 
        k_th_phi) = gen_ks(pname,hname,equatorial)
    for r in ratios:
      if freq:
        for x in bbh.invert(k_r_phi,om_phi,r):
          r_phi_crossings.append(x)
        for x in bbh.invert(k_r_th,om_phi_trimmed,r):
          r_th_crossings.append(x)
        for x in bbh.invert(k_th_phi,om_phi_trimmed,r):
          th_phi_crossings.append(x)
      else:
        for x in bbh.invert(k_r_phi,om_r_t,r):
          r_phi_crossings.append(x)
        for x in bbh.invert(k_r_th,trim_time,r):
          r_th_crossings.append(x)
        for x in bbh.invert(k_th_phi,trim_time,r):
          th_phi_crossings.append(x)
    return r_phi_crossings,r_th_crossings,th_phi_crossings 


def ks(pname="Omega.pkl",hname="Horizons.h5",equatorial=False,freq=False):
  import cPickle as pickle
  import bbh_processing.utilities as bbh
  import numpy as np
 
  ratios = [1., 0.5, 2./3., 3./4., 4./5.] 
  ratlabs = ["1:1", "1:2", "2:3", "3:4", "4:5"]
  for z in zip(ratios,ratlabs):
    plt.axhline(y=z[0],label=z[1])
  if equatorial:
    om_phi, om_r_t, k_r_phi = gen_ks(pname,hname,equatorial) 
    #no trimming necessary
    out=file("ks_equatorial.dat",'w')
    out.write("""#
      # [1] = t
      # [2] = om_r/om_phi 
      """)
    np.savetxt(f,np.transpose([om_r_t,k_r_phi])) 
    if freq:
      the_plot = plt.plot(om_phi,k_r_phi) 
      plt.xlabel("Om_phi")
    else:
      the_plot = plt.plot(om_r_t,k_r_phi)
      plt.xlabel("Time(M)")
    
    lab=['r:phi']
    plt.legend(lab, loc='lower left',fontsize=8)
    plt.ylabel("Ratio")
    plt.title("Frequency Ratios")
    return the_plot

  else:
    (om_phi,om_phi_trimmed,om_r_t, k_r_phi, trim_time, k_r_th, trim_time, 
        k_th_phi) = gen_ks(pname,hname,equatorial)
    
    if freq:
      the_plot = plt.plot(om_phi,k_r_phi,
                        om_phi_trimmed,k_r_th,
                        om_phi_trimmed,k_th_phi)
      plt.xlabel("Om_phi")
    else:
      the_plot = plt.plot(om_r_t,k_r_phi,
                        trim_time,k_r_th,
                        trim_time,k_th_phi)
      plt.xlabel("Time(M)")
    lab=['r:phi','r:th','th:phi']
    plt.legend(lab, loc='lower left',fontsize=8)
    plt.ylabel("Ratio")
    plt.title("Frequency Ratios")
    return the_plot



#Plot E and E_dot
def energy(fname="Fluxes.dat"):
  import bbh_processing.ffile_tools as ff 
  T = ff.t(fname)
  E = ff.e(fname)
  Edot = ff.de_dt(fname)
  the_plot = plt.plot(T,E)
  plt.xlabel("Time(M)")
  plt.ylabel("E")
  plt.title("Energy")
  return the_plot

def edot(fname="Fluxes.dat"):
  import bbh_processing.ffile_tools as ff 
  T = ff.t(fname)
  Edot = ff.de_dt(fname)
  
  the_plot = plt.plot(T,Edot)
  plt.xlabel("Time(M)")
  plt.ylabel("Edot")
  plt.title("Flux")
  return the_plot

#Plot s as a function of time
def s_t(hname="Horizons.h5",equatorial=False):
  import bbh_processing.utilities as bbh
  import bbh_processing.hfile_tools as hf 
  t = hf.t(hname)
  S = hf.chi_inertial(hname)
  S_mag = la.norm(S,axis=1)
  
 
  lab=['S_x','S_y','S_z','|S|']
  if equatorial:
    lab=['S_x','S_y','|S|']
  the_plot = plt.plot(t,S,t,S_mag)
  plt.legend(lab, loc='upper left',fontsize=8)
  plt.xlabel("Time(M)")
  plt.ylabel("S")
  plt.title("ChiInertial")
  return the_plot

#Plot J projected onto S. This involves an interpolation.
def j_on_s(hname="Horizons.h5",fname="Fluxes.dat",equatorial=False):
  import bbh_processing.utilities as bbh
  import bbh_processing.ffile_tools as ff 
  import bbh_processing.hfile_tools as hf 
  t = hf.t(hname)
  T = ff.t(fname)
  S = hf.chi_inertial(hname)
  J = bbh.interp_onto(ff.j(fname),T,t,'slinear')
  J_S = bbh.projected_onto(J,S)
   
  J_S_mag = la.norm(J_S,axis=1) 
 
  lab=['x','y','z','mag']
  if equatorial:
    lab=['x','y','mag']
  the_plot = plt.plot(t,J_S,t,J_S_mag)
  plt.legend(lab,loc='upper left',fontsize=8)
  plt.xlabel("Time(M)")
  plt.ylabel("J")
  plt.title("J projected onto S")
  return the_plot

#Plot J. No interpolation.
def j(fname="Fluxes.dat"):
  import bbh_processing.utilities as bbh
  import bbh_processing.ffile_tools as ff 
  
  T = ff.t(fname)
  J = ff.j(fname)
  J_mag = la.norm(J,axis=1)

  lab=['J_x','J_y','J_z','|J|']
  the_plot = plt.plot(T,J,T,J_mag)
  plt.legend(lab,loc='upper left',fontsize=8)
  plt.xlabel("Time(M)")
  plt.ylabel("J")
  plt.title("J in simulation coordinates")
  return the_plot

#Plot P projected onto S. This involves an interpolation.
def p_on_s(hname="Horizons.h5",fname="Fluxes.dat",equatorial=False):
  import bbh_processing.utilities as bbh
  import bbh_processing.ffile_tools as ff 
  import bbh_processing.hfile_tools as hf 
  t = hf.t(hname)
  S = hf.chi_inertial(hname)
  T = ff.t(fname)
  P = bbh.interp_onto(ff.p(fname),T,t,'slinear')
  P_S = bbh.projected_onto(P,S)
  P_S_mag = la.norm(P_S, axis=1) 
  
  lab=['x','y','z','mag']
  if equatorial:
    lab=['x','y','mag']
  the_plot = plt.plot(t,P_S,t,P_S_mag)
  plt.legend(lab,loc='upper left',fontsize=8)
  plt.xlabel("Time(M)")
  plt.ylabel("P")
  plt.title("P projected onto S")
  return the_plot

#Plot P. No interpolation.
def p(fname="Fluxes.dat"):
  import bbh_processing.utilities as bbh
  import bbh_processing.ffile_tools as ff 
  T = ff.t(fname)
  P = ff.p(fname)
 
  P_mag = la.norm(P,axis=1) 
  
  lab=['P_x','P_y','P_z','|P|']
  the_plot = plt.plot(T,P,T,P_mag)
  plt.legend(lab,loc='upper left',fontsize=8)
  plt.xlabel("Time(M)")
  plt.ylabel("P")
  plt.title("P in simulation coordinates")
  return the_plot

#Plot theta.
def theta_t(hname="Horizons.h5"):
  import bbh_processing.hfile_tools as hf
  T = hf.t(hname)
  theta = hf.theta_r_chi(hname)

  the_plot = plt.plot(T,theta)
  plt.xlabel("Time(M)")
  plt.ylabel("Radians")
  plt.title("Theta (r and chi)")

#Plot the components of rhat projected onto chi.
def rhat_projected(hname="Horizons.h5"):
  import bbh_processing.hfile_tools as hf
  import bbh_processing.utilities as bbh
  T = hf.t(hname)
  chi_hat = hf.chi_hat(hname)
  r_vec = hf.r(hname)
  r_proj = bbh.projected_onto(r_vec,chi_hat)
  r_proj_hat = bbh.unit_vector_time_series(r_proj)
  lab=['x','y','z']
  the_plot = plt.plot(T,r_proj_hat)
  plt.xlabel("Time(M)")
  plt.ylabel("Length(M)")
  plt.title("Separation vector projected onto spin axis")
  return the_plot
#plt.show()
