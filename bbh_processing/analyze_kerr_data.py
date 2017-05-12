"""
This bundles methods to process data from Aaron's Mathematica notebook.
"""
import numpy as np
import numpy.linalg as linalg
from numpy.lib import stride_tricks
#import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import utilities as bbh
from scipy import optimize
import bbh_processing.unwrap as unwrap
from scipy.signal import savgol_filter
import peak_finders as peaks

# def interpolated_grid_coords(x_raw):
    # """
    # Take raw data from the .dat file and return data interpolated to an even
    # timestep.
    # """
    # time = np.copy(x_raw[:, 0])
    # x_interp = x_raw[:, 1]
    # y_interp = x_raw[:, 2]
    # z_interp = x_raw[:, 3]
    # xspl = UnivariateSpline(time, x_interp, k=4, s=0)
    # yspl = UnivariateSpline(time, y_interp, k=4, s=0)
    # zspl = UnivariateSpline(time, z_interp, k=4, s=0)
    # tnew = np.arange(time[0], time[-1], 0.5)
    # #tnew = np.arange(t[0], t[-1], (t[1]-t[0])/10.)
    # xnew = xspl(tnew)
    # ynew = yspl(tnew)
    # znew = zspl(tnew)
    # return np.column_stack((tnew, xnew, ynew, znew))

def rmag(rvec):
    """
    Compute the magnitude of r
    """
    #time = np.copy(rvec[:, 0])
    rdat = np.sqrt(rvec[:, 0]**2 + rvec[:, 1]**2 + rvec[:, 2]**2)
    #rdat = linalg.norm((rvec[:, 0], rvec[:, 1], rvec[:, 2]), ord=2, axis=0)
    return rdat
    #return np.column_stack((time, rdat))

def rdot(time, rvec):
    """
    Return the gradient of r.
    """
    dtime = np.gradient(time)
    dxdt = np.gradient(rvec[:, 0], dtime)
    dydt = np.gradient(rvec[:, 1], dtime)
    dzdt = np.gradient(rvec[:, 2], dtime)
    return np.column_stack((dxdt, dydt, dzdt))

def theta(rvec, chi):
    """
    Return the angle between chi and rvec.
    """
    chi_ref = stride_tricks.as_strided(np.array(chi), strides=(0, 1*8), 
                                       shape=rvec.shape)
    return bbh.angle(rvec, chi_ref)

def costheta(rvec, chi):
    """
    Return the angle between chi and rvec.
    """
    chi_ref = stride_tricks.as_strided(np.array(chi), strides=(0, 1*8), 
                                       shape=rvec.shape)
    return bbh.cosang(rvec, chi_ref)
 
def rho(rvec):
    return np.sqrt(rvec[:, 0]**2 + rvec[:, 1])**2

def rhoref(rvec):
    shp = (len(rvec), 2)
    rhovec = np.zeros(shp)
    rhovec[:, 0] = rvec[:, 0]
    rhovec[:, 1] = rvec[:, 1]
    refval = rhovec[0, :]
    ref = np.tile(refval, (len(rvec), 1)) 
    return rhovec, ref

def phi(rvec):
    rhovec, ref = rhoref(rvec)
    return bbh.angle(rhovec, ref)     

def cosphi(rvec):
    rhovec, ref = rhoref(rvec)
    return bbh.cosang(rhovec, ref)     

def unwrapphi(rvec):
    thisphi = phi(rvec)
    unwrapped, flipped = unwrap.unwrap(thisphi, junk_inds=0)
    return np.unwrap(flipped)

def loadfile(fname="KerrAnalyticOrbitData.dat"):
    """
    Read the file.
    """
    with open(fname, 'r') as myfile:
        xraw = np.loadtxt(myfile)
    return xraw

def even_sampled_omega_npts(orbpath, npts, method="savgol",
                       return_cos=False):
    rraw = loadfile(orbpath)
    time = np.linspace(rraw[0, 0], rraw[-1, 0], npts)
    x = UnivariateSpline(rraw[:, 0], rraw[:, 1], k=4, s=0)(time)
    y = UnivariateSpline(rraw[:, 0], rraw[:, 2], k=4, s=0)(time)
    z = UnivariateSpline(rraw[:, 0], rraw[:, 3], k=4, s=0)(time)
    rvec = np.column_stack([x, y, z])
    
    ommag = totalfreq(time, rvec, method=method)
    if not return_cos:
        return time, ommag
    else:
        costhetaraw = costheta(rraw[:, 1:], [0, 0, 1])
        costh = UnivariateSpline(rraw[:, 0], costhetaraw, k=4, s=0)(time)
        rplane = np.copy(rvec)
        rplane[:, 2] = 0.
        omplane = savgol_omegafromr(time, rplane, 7, delta=(time[1]-time[0]))
        return time, ommag, costh, omplane
      
def even_sampled_omega(orbpath, ptsperorbit=1000, method="savgol",
                       return_cos=False):
    rraw = loadfile(orbpath)
    rmagraw = np.sqrt(rraw[:, 1]**2 + rraw[:, 2]**2 + rraw[:, 3]**2)
    pt, pv = peaks.conservative_peaks(rraw[:, 0], rmagraw, refine=False)
    norbits = len(pt)
    print "Norbits: ", norbits
    npts = norbits * ptsperorbit
    return even_sampled_omega_npts(orbpath, npts, method=method, return_cos=return_cos)
    
def periodogramfft(d, time, data):
    from scipy.signal import hann, periodogram
    import numpy.fft as fft
    print "d=", d
    
     
    pidx, pv, tidx, tv = peaks.turning_points(data, minidx=0)
    tts = time[pidx]
    rspl = UnivariateSpline(time, data, k=4, s=0)
    t0 = tts[1]
    tf = tts[d+1]
    dt = time[1]-time[0]

    
    twin = np.arange(t0, tf, dt)
    N = len(twin) 
    fs = 1/dt
    #twin1 = np.arange(t0, tts[3], dt)
    #Nfft = len(twin1)
    #fs = N/(tf-t0)
    Nfft=100000
    rwin = rspl(twin)
    rwin = rwin - np.mean(rwin)
    hwin = hann(len(twin))
    rwin = rwin*hwin
    #f, pgram = periodogram(hwin, fs=fs, nfft=Nfft)
    fftout = fft.rfft(rwin, 50*len(rwin))
    fftfreq = fft.rfftfreq(50*len(rwin), dt)
    print len(fftout)
    print len(fftfreq)
    return 2*np.pi*fftfreq, np.abs(fftout)
  
def periodogram(d, time, data):
    from scipy.signal import hann, lombscargle
    print "d=", d
    #datamean = np.mean(data)
    datashift = data - np.mean(data) 
    pidx, pv, tidx, tv = peaks.turning_points(datashift, minidx=0)
    tts = time[pidx]
    rspl = UnivariateSpline(time, datashift, k=4, s=0)
    twin = np.linspace(tts[1], tts[d+1], 5000)
    rwin = rspl(twin)
    hwin = hann(len(twin))
    rwin = rwin*hwin
    f = np.linspace(0.005, 0.1, 3000)
    pgram = lombscargle(twin, rwin, f)
    pgram = np.sqrt(4*(pgram/rwin.shape))
    return f, pgram

def savgol_omegafromr(t, r, N, polorder=3, delta=0.5):
    r2 = r[:, 0]**2 + r[:, 1]**2 + r[:, 2]**2
    dv = savgol_filter(r, N, polorder, deriv=1, delta=delta, axis=0)  
    omvals = np.cross(r, dv)/r2[:, None]
    ommag = np.sqrt(omvals[:, 0]**2 + omvals[:, 1]**2 + omvals[:, 2]**2) 
    #OmegaMag = np.column_stack((t, ommag))
    return ommag

def totalfreq(time, r_in, method="savgol"):
    if method=="savgol":
        return savgol_omegafromr(time, r_in, 7, delta=(time[1] - time[0]))
    elif method=="direct":
        return totalfreqdirect(time, r_in)
    elif method=="spline":
        return totalfreqspline(time, r_in)

def totalfreqspline(time, r_in):
    rsqr = (rmag(r_in))**2
    xspl = UnivariateSpline(time, r_in[:, 0], k=4, s=0)
    yspl = UnivariateSpline(time, r_in[:, 1], k=4, s=0)
    zspl = UnivariateSpline(time, r_in[:, 2], k=4, s=0)

    xdot = xspl.derivative()(time)
    ydot = yspl.derivative()(time)
    zdot = zspl.derivative()(time)
    r_dot = np.column_stack((xdot, ydot, zdot))

    cross = np.cross(r_in, r_dot)
    quot = np.zeros(cross.shape)
    for i in range(0, 3):
        quot[:, i] = np.divide(cross[:, i], rsqr)
    mag = np.sqrt(quot[:, 0]**2 + quot[:, 1]**2 + quot[:, 2]**2)
    return mag

def totalfreqdirect(time, r_in):
    """
    Compute the orbital frequency vector.
    """
    #time = np.copy(r_in[:, 0])
    rsqr = (rmag(r_in))**2
    r_dot = rdot(time, r_in)
    #print r_dot.shape
    cross = np.cross(r_in, r_dot)
    quot = np.zeros(cross.shape)
    for i in range(0, 3):
        quot[:, i] = np.divide(cross[:, i], rsqr)
    mag = np.sqrt(quot[:, 0]**2 + quot[:, 1]**2 + quot[:, 2]**2)
    return mag


def polyderiv(time, data, N=5):
    times = []
    ranges = []
    length = len(data)
    if length != len(time):
        raise ValueError("Inconsistent sizes.")

    for i in range(N, length-N):
        l = i-N
        r = i+N+1
        times.append(time[l:r])
        ranges.append(data[l:r])
    x = []
    v = []
    for t, r in zip(times, ranges):
        poly = np.polyfit(t, r, 2)
        x.append(poly[0])
        v.append(poly[1])
     
    return np.array(x), np.array(v)

def fitfreq(time, r_in):
    """
    Compute the orbital frequency vector.
    """
    #time = np.copy(r_in[:, 0])
    xfit = []
    vfit = []
    for i in range(0, 3):
        xf, vf = polyderiv(time, r_in[:, i])
        xfit.append(xf)
        vfit.append(vf)
    xarr = np.transpose(np.array(xfit))
    varr = np.transpose(np.array(vfit))
    r_mag = rmag(xarr)
    rsqr = r_mag**2
    omega = np.cross(xarr, varr)
    for i in range(0, 3):
        omega[:, i] = omega[:, i] / rsqr
    omegamag = np.copy(omega[:, 0:2])
    omegamag[:, 1] = np.sqrt(omega[:, 0]**2 + omega[:, 1]**2 + omega[:, 2]**2)
    return omegamag
    
def makeabdulfreq(rraw, N=10):
    """Given numpy arrays for the black hole locations xA, xB,
    perform fits to 2N+1 data-points around each point, and from
    the fit compute the instantaneous orbital frequency,
            Omega = r\times \dot{r} / r^2
    return Omega, OmegaMag
"""
    # compute orbital frequency as a function of time
    fitfunc = lambda p,x: p[0]+p[1]*x+p[2]*x**2
    errfunc = lambda p,x,y: fitfunc(p,x)-y

    def FitData(data, N):
        """given a numpy array data with time as first column, perform fits 
        covering N points before and after each data-point for each column.
        return the fitted values and their first time-derivatives as a 
        numpy arrays, with first column time"""
        # collect output data, first column being time
        last_idx=len(data)-N-1
        x_fit=[ list(data[N:last_idx,0]) ]
        v_fit=[ list(data[N:last_idx,0]) ]
        for ycol in range(1,data.shape[1]):
            p0=[data[N,ycol],0,0]
            x_tmp=[]
            v_tmp=[]
            for idx in range(N, len(data)-N-1):
                T=data[idx,0]
                x=data[idx-N:idx+N+1,0]-T # shift target to zero
                y=data[idx-N:idx+N+1,ycol]
                p0,success=optimize.leastsq(errfunc, p0[:],args=(x,y))
                x_tmp.append(p0[0])
                v_tmp.append(p0[1])
            x_fit.append(x_tmp)
            v_fit.append(v_tmp)
        return np.transpose(np.array(x_fit)),\
               np.transpose(np.array(v_fit))

    r_fit ,v = FitData(rraw, N)

    #compute Orbital frequency
    #dr=xA_fit-xB_fit
    #dv=vA-vB
    r2=(r_fit[:,1]**2+r_fit[:,2]**2+r_fit[:,3]**2)
    Omega = r_fit
    Omega[:, 1]=(r_fit[:,2]*v[:,3]-r_fit[:,3]*v[:,2])/r2
    Omega[:, 2]=(r_fit[:,3]*v[:,1]-r_fit[:,1]*v[:,3])/r2
    Omega[:, 3]=(r_fit[:,1]*v[:,2]-r_fit[:,2]*v[:,1])/r2
    OmegaMag = Omega[:, 0:2].copy()
    OmegaMag[:, 1]=np.sqrt(Omega[:, 1]**2+Omega[:, 2]**2+Omega[:, 3]**2)
    return Omega, OmegaMag

if __name__ == "__main__":
    raise NotImplementedError("This module cannot be run from command line.")
