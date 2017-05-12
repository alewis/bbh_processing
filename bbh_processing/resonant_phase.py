import numpy as np
from scipy import integrate
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import peak_finders as pk


def time_of_phase(time, phases, t_res, p_res, diff):
    shift = np.abs(phases - p_res) - diff
    spl = UnivariateSpline(time, shift, k=3, s=0)
    rts = spl.roots()
    if len(rts)<2:
      return -1, -1
    else:
      arg = np.argsort(np.abs(rts-t_res))
      t1 = rts[arg[0]]
      t2 = rts[arg[1]]
      if t1 >= t2:
        t_l = t2
        t_r = t1
      else:
        t_l = t1
        t_r = t2
      return t_l, t_r

def resonance_plotdata(t, omph, om_small, om_big, k_lesser, k_bigger, 
                       tr=None, r=None, tth=None, th=None,
                       Npts=10000, resl=1):
    integrand = k_bigger*om_small - k_lesser*om_big
    tint = np.linspace(t[0], t[-1], Npts) 
    omphint = np.linspace(omph[0], omph[-1], Npts)
    intspl = UnivariateSpline(t, integrand, k=3, s=0)
    t_res = np.mean(intspl.roots())

    integrand_int = intspl(tint) 
    phase = integrate.cumtrapz(integrand_int, x=tint, initial=0.)
    phase_res = UnivariateSpline(tint, phase, k=4, s=0)(t_res)
    t_left, t_right = time_of_phase(tint, phase, t_res, phase_res, resl)
  
    cosphase = np.cos(phase)        
    cosspl = UnivariateSpline(tint, cosphase, k=4, s=0)
    at_points = cosspl(t)

    omphspl = UnivariateSpline(tint, omphint, k=1, s=0)
    omph_left = omphspl(t_left)
    omph_res = omphspl(t_res)
    omph_right = omphspl(t_right)

    if tr is not None:
        Nr = orbits_on_resonance(tr, r, t_left, t_right)
        Nth = orbits_on_resonance(tth, th, t_left, t_right, refine=False) 
        traj_unref = resonant_trajectory(tr, r, th, t_left, t_right, refinex1=False)
        traj_ref = resonant_trajectory(tr, r, th, t_left, t_right, refinex1=True)
    else:
        Nr = -1
        Nth = -1
        traj_unref = np.array([-1, -1])
        traj_ref = np.array([-1, -1])
    
    return [(t_left, t_res, t_right), 
            (omph_left, omph_res, omph_right),
            (tint, omphint, cosphase), 
            (t, omph, at_points),
            (Nr, Nth),
            traj_unref,
            traj_ref]

def orbits_on_resonance(t, r, t_left, t_right, refine=True):
    pt, pv = pk.conservative_peaks(t, r, refine=refine, minidx=200)
    argl = np.argmin(np.abs(pt - t_left))   
    argr = np.argmin(np.abs(pt - t_right))
    return argr-argl
    
def resonant_trajectory(t, x1, x2, t_left, t_right, Npts=10000, refinex1=False,
                        t_res=None): 
    if refinex1:
        pt, pv, p_ref, pv_unref = pk.conservative_peaks(t, x1, returnrefined=True) 
        x1spl = UnivariateSpline(t, p_ref, k=4, s=0)
    else:
        x1spl = UnivariateSpline(t, x1, k=4, s=0)
    x2spl = UnivariateSpline(t, x2, k=4, s=0)
    t_traj = np.linspace(t_left, t_right, Npts)
    x1_traj = x1spl(t_traj)
    x2_traj = x2spl(t_traj)
    return [x1_traj, x2_traj]

def trajpassage(N, traj):
    """Compute from the trajectory the index of the resonance (i.e. the halfway
    index) and those of N full passages in either direction from it.
    """
    residx = len(traj)/2
    pidx, pv = pk.turning_points(traj, joined=True)
    pidx = np.array(pidx, dtype=int)
    try:
        right_ex = np.argmax(pidx>residx)
    except ValueError:
        return [-1, residx, -1]
    right_done = right_ex + N
    if right_done >= len(pidx):
        right_idx = len(traj)-1
    else:
        right_idx = pidx[right_done]
    
    left_ex = right_ex - 1
    left_done = left_ex - N
    if left_done < 0:
        left_idx = 0
    else:
        left_idx = pidx[left_done]
    # plt.plot(traj)
    # plt.plot(residx, traj[residx], marker="o", markersize=10)
    # plt.plot(left_idx, traj[left_idx], marker="o", markersize=10)
    # plt.plot(right_idx, traj[right_idx], marker="o", markersize=10)
    # plt.show()
    return [left_idx, residx, right_idx]

def toruspassage(N1, traj1, N2, traj2):
    left1, res1, right1 = trajpassage(N1, traj1)
    left2, res2, right2 = trajpassage(N2, traj2)

    left = left1
    if left2 < left:
        left = left2

    right = right1
    if right2 > right:
        right = right2
    
    return [left, res1, right]



#Compute 
#P=INT[n*Omega_i + k*Omega_j]dt where the n and k specify the resonant order. This is used in Hughes 2014 to diagnose passage through a resonance.
def resonant_phase(t, om_small, om_big, ksmall, kbig):

    integrand = kbig*om_small - ksmall*om_big
    phase = integrate.cumtrapz(integrand, x=t, initial=0.)
    return phase 

#Just the cosine of the resonant phase.
def cosine_phase(t, om_i, om_j, n, k):
    phase = resonant_phase(t, om_i, om_j, n, k)
    return np.cos(phase)

def tres(t, om_small, om_big, n, k):
    res = float(n)/float(k)
    integrand = k*om_small - n*om_big
    kij_int = UnivariateSpline(t, integrand, k=3, s=0)
    rts = kij_int.roots()
    return np.mean(rts)
    #theroot = np.argmin(np.abs(rts-tguess))

#compute the resonant time, when 'phase' (not the cosine phase) attains
#its maximum.
def t_resonance(phase):
  turns, turn_vals = pk.turning_points(phase[:,1])
  i = np.argmax(turn_vals)
  return phase[turns[i],0]
  #t_res=14010
  #return t_res

# #compute the time spent on resonance as 2/sqrt(k*dq1/dt + n*dq2/dt)
# def time_on_resonance(q1,q2,k1,k2,phase,t_res):
  # k1_q1dot_spl = interp.UnivariateSpline(q1[:,0],k1*q1[:,1],k=4,s=0).derivative()
  # k2_q2dot_spl = interp.UnivariateSpline(q2[:,0],k2*q2[:,1],k=4,s=0).derivative()
  # return 2./np.sqrt(np.fabs(k1_q1dot_spl(t_res) + k2_q2dot_spl(t_res)))

# #compute the time spent on resonance as the time when the phase changes from its
# #resonance value by no more than one radian
# def time_on_resonance_xrad(phase,t_res,x=1,tstep=0.5):
  # phase_spl = interp.UnivariateSpline(phase[:,0],np.arccos(phase[:,1]),k=3,s=0) 
  # p_res = phase_spl(t_res)
  # #the left value 
  # ti=t_res
  # while True:
   # ti -= tstep
   # if np.fabs(phase_spl(ti) - phase_spl(t_res)) > x:
     # break

  # #the right value 
  # tf=t_res
  # while True:
   # tf += tstep 
   # if np.fabs(phase_spl(tf) - phase_spl(t_res)) > x:
     # break
  # return ti,tf 

# def trajectory_around_resonance_xrad(q1,q2,phase,t_res,x=0.5):
  # ti,tf = time_on_resonance_xrad(phase,t_res,x=0.5)
  # q1spl = interp.UnivariateSpline(q1[:,0],q1[:,1],k=3,s=0)
  # q2spl = interp.UnivariateSpline(q2[:,0],q2[:,1],k=3,s=0)
  
  # #column 0 - time
  # #column 1 - q1 
  # #column 2 - q2
  # sz = 2.*(tf-ti)
  # traj = np.zeros((sz,3)) 
  # trng = np.linspace(ti,tf,sz) 
  # traj[:,0]=np.copy(trng)
  # traj[:,1]=q1spl(trng)
  # traj[:,2]=q2spl(trng)
  # return traj
# def trajectory_around_resonance(q1,q2,k1,k2,phase,t_res,N=50):
  # T = time_on_resonance(q1,q2,k1,k2,phase,t_res)
  # q1spl = interp.UnivariateSpline(q1[:,0],q1[:,1],k=3,s=0)
  # q2spl = interp.UnivariateSpline(q2[:,0],q2[:,1],k=3,s=0)
  
  # #column 0 - time
  # #column 1 - q1 
  # #column 2 - q2
  # sz = 2.*T
  # traj = np.zeros((sz,3)) 
  # trng = np.linspace(t_res-T/2.,t_res+T/2.,sz) 
  # traj[:,0]=np.copy(trng)
  # traj[:,1]=q1spl(trng)
  # traj[:,2]=q2spl(trng)
  # return traj
