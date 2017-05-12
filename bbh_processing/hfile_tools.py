import h5py
import numpy as np
import utilities as bbh
import math
import numpy.linalg as la
import cartesian_basis_series as basis
import peak_finders as peaks
from numpy.lib import stride_tricks
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter
#Functions to process things from Horizons.h5 files.
###########################################################
#Functions to directly extract things.
###########################################################

#Convenience function to unify coordinate retrieval.
def get_coord(name, path, rotate=False):
    """Retrieves the named coord. If rotate is true, rotates into the frame
       in which chi is parallel with the z axis.
    """
    r_array = r(path, rotate)
    #if rotate==False:
    chi_array = chi_inertial(path, rotate) 
    # else:
        # zhat = np.array([0., 0., 1.])
        # chi_array = stride_tricks.as_strided(zhat,strides=(0,1*8),
                      # shape=r_array.shape)
    if name == "t":
        coord = t(path)
    elif name == "xA":
        f=h5py.File(path, 'r')
        return f["/AhA.dir"]["CoordCenterInertial.dat"].value[:,1:]
    elif name == "xB":
        f=h5py.File(path, 'r')
        return f["/AhB.dir"]["CoordCenterInertial.dat"].value[:,1:]
    elif name == "rvec":
        coord = r_array
    elif name == "x":
        coord = r_array[:, 0]
    elif name == "y":
        coord = r_array[:, 1]
    elif name == "z":
        coord = r_array[:, 2]
    elif name == "r":
        coord = la.norm(r_array, axis=1)
    elif name == "r2":
        coord = r_array[:, 0]**2 + r_array[:, 1]**2 + r_array[:, 2]**2 
    elif name == "rho":
        rho = r_array
        rho[:, 2] = 0.
        coord = rho
    elif name == "chi":
        coord = chi_array
    elif name == "chimag":
        coord = la.norm(chi_array, axis=1)
    elif name == "chihat":
        coord = np.zeros(chi_array.shape)
        chimag = la.norm(chi_array, axis=1)
        coord[:, 0] = chi_array[:, 0] / chimag
        coord[:, 1] = chi_array[:, 1] / chimag
        coord[:, 2] = chi_array[:, 2] / chimag
    elif name == "phi":
        coord = np.arccos(cosphi(r_array))
    elif name == "cosphi":
        coord = cosphi(r_array)
    elif name == "phi_phase":
        coord = bbh.accumulating_phi(r_array) 
    elif name == "costheta":
        coord = bbh.cosang(chi_array, r_array)
    elif name == "theta":
        coord = np.arccos(bbh.cosang(chi_array, r_array))
    elif name == "theta_phase":
        coord = bbh.phase(r_array, chi_array)
    else:
        raise ValueError("Invalid coordinate name '"+name+"'")
    return coord


def savgol_omega(path, N=7, polorder=3, delta=0.5):
    t = get_coord("t", path)
    xA = np.column_stack((t, get_coord("xA", path)))
    xB = np.column_stack((t, get_coord("xB", path)))
    Omega, OmegaMag = savgol_omegafunc(xA, xB, N, polorder=polorder, delta=delta)
    peaklist, smooth = peaks.turning_points(OmegaMag[:, 1], return_smooth=True)
    
    OmegaMag[:, 1] = smooth
    #smoothed = smooth_savgol(OmegaMag)
    #OmegaMag_smooth = smooth_omegamag(OmegaMag)
    return Omega, OmegaMag 

#def smooth_omegamag(data):

def rperp(r, chi):
    """r_perp = r - hat(s) (hat(s) . r)
    """
    chihat = bbh.unit_vector(chi)
    chidotr = bbh.timedot(chihat, r)
    chiproj = bbh.scaled_vector(chidotr, chihat)
    rperp = r - chiproj
    return rperp

def savgol_omegafunc(xA, xB, N, polorder=3, delta=0.5):
    dx = xA[:, 1:] - xB[:, 1:]
    return savgol_omegafromr(xA[:, 0], dx, N, polorder=polorder, delta=delta)

def savgol_omegafromr(t, r, N, polorder=3, delta=0.5):
    r2 = r[:, 0]**2 + r[:, 1]**2 + r[:, 2]**2
    dv = savgol_filter(r, N, polorder, deriv=1, delta=delta, axis=0)  
    omvals = np.cross(r, dv)/r2[:, None]
    Omega = np.column_stack((t.copy(), omvals))
    
    ommag = np.sqrt(omvals[:, 0]**2 + omvals[:, 1]**2 + omvals[:, 2]**2) 
    OmegaMag = np.column_stack((t, ommag))
    return Omega, OmegaMag
    

#Angle between rho and its initial value.
#Rho should be the x and y components of r after rotation into the chi frame.
def cosphi(r_array):
    rho = r_array[:, 0:2]
    rho_ref = stride_tricks.as_strided(rho[0,:],strides=(0,1*8),shape=rho.shape)
    return bbh.cosang(rho,rho_ref)

#Return an array of the time part of the coord-centre time series
def t(fname="Horizons.h5"):
  f=h5py.File(fname,'r')
  xA = np.array(f["/AhA.dir"]["CoordCenterInertial.dat"])
  f.close()
  return xA[:,0] 

#Return the spin vector
def chi_inertial(fname="Horizons.h5", rotate=False):
  #main
  f=h5py.File(fname,'r')
  #spinA = f["/AhA.dir"]["DimensionfulInertialSpin.dat"].value[:,1:]
  spinA = f["/AhA.dir"]["chiInertial.dat"].value[:,1:]
  if rotate:
      t = get_coord("t", fname)
      chi = chi_inertial(fname)
      chibasis = basis.CartesianBasisSeries(t, chi)
      spinA = chibasis.rotate(spinA)
  f.close()
  return spinA 


#Return the separation vector (r)
def r(fname="Horizons.h5", rotate=False):
  #main
  f=h5py.File(fname,'r')
  xA = f["/AhA.dir"]["CoordCenterInertial.dat"].value
  xB = f["/AhB.dir"]["CoordCenterInertial.dat"].value
  f.close()
  dr=xA[:, 1:]-xB[:, 1:]
  
  if rotate:
      #t = get_coord("t", fname)
      t = xA[:, 0]
      chi = chi_inertial(fname)
      chibasis = basis.CartesianBasisSeries(t, chi)
      dr = chibasis.rotate(dr)
      #print r
  return dr 

#Return the angle between the spin and separation vector
def costheta_r_chi(r_array,  chi_array):
  #main
  spinA = chi_inertial(fname)
  return bbh.cosang(spinA, r_array)

#Unit separation vector
def r_hat(fname="Horizons.h5"):
  return bbh.unit_vector(r(fname))

#Unit spin vector
def chi_hat(fname="Horizons.h5"):
  return bbh.unit_vector(chi_inertial(fname))

#The xy-position of the centre of mass as a function of time.
def axial_cm_motion(q,fname="Horizons.h5",mtot=1):
  f=h5py.File(fname,'r')
  m1 = mtot*(1-1/(q+1))
  m2 = mtot/(q+q)
  xA = np.array(f["/AhA.dir"]["CoordCenterInertial.dat"])
  xB = np.array(f["/AhB.dir"]["CoordCenterInertial.dat"])
  #coordinate components
  Xa = m1*xA[:,1]/mtot
  Ya = m1*xA[:,2]/mtot
  Xb = m2*xB[:,1]/mtot
  Yb = m2*xB[:,2]/mtot
  
  xcm = Xa-Xb
  ycm = Ya-Yb
  return (xcm, ycm)
  #Xa -= xcm[0]
  #Xa -= ycm[0]
  #Xa -= zcm[0]

#The eccentricity from e= \frac{r_a - r_p}{r_a + r_p}
def ecc_from_app_peri(newtimes=None, fname="Horizons.h5"):
    rmag = get_coord("r", fname)
    t = get_coord("t", fname)
    t_a, r_a, t_p, r_p = peaks.peaks_from_polyfit(t, rmag)
    ra_spl = UnivariateSpline(t_a, r_a, s=0)
    rp_spl = UnivariateSpline(t_p, r_p, s=0)

    if newtimes is not None:
        inttimes = newtimes
    else:
        inttimes = [(ta+tp)/2. for ta, tp in zip(t_a, t_p)]
    ra_int = ra_spl(inttimes)
    rp_int = rp_spl(inttimes)
    ecc = (ra_int - rp_int) / (ra_int + rp_int)
    # ecclam = lambda ra, rp: (ra - rp)/(r_a + r_p)
    # ecc = [ecclam(ra, rp) for ra, rp in zip(ra_mid, rp_mid)] 


    return inttimes, ecc





