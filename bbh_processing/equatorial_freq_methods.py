from scipy.interpolate import UnivariateSpline
from scipy.integrate import simps
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import cPickle as pickle
import hfile_tools as hfile
import peak_finders as peak
##############################################################################
#INTERFACE
##############################################################################
def get_peaks(time, data, pktype="peak",refine=True, refine_pktype="peak", Npk=5,
              polyorder=3, splorder=3, minidx=200, peak_method="spline"):
    if pktype == "joined":
        pt, pv = peak.joined_peaks(time, data, refine=refine,
                                   refine_pktype=refine_pktype, Npk=Npk,
                                   polyorder=polyorder, splorder=splorder,
                                   minidx=minidx, peak_method=peak_method)
        joined=True
    else:
        pt, pv = peak.conservative_peaks(time, data, refine=refine,
                                         refine_pktype=refine_pktype, Npk=Npk,
                                         polyorder=polyorder, splorder=splorder,
                                         minidx=minidx, peak_method=peak_method)
        joined=False
    return pt, pv, joined


def ecc_sqrt(pt, pv, tt, tv, int_t, minidx=200):
    pspl = UnivariateSpline(pt, pv, k=5, s=0)
    tspl = UnivariateSpline(tt, tv, k=5, s=0)
     
    r_p = np.sqrt(pspl(int_t))
    r_a = np.sqrt(tspl(int_t))
    ecc = (r_p - r_a)/(r_p + r_a)
    return np.array(ecc)

# def ecc_sqrt(omtime, ommag, int_t, minidx=200):
    # pt, pv = peak.conservative_peaks(omtime, ommag, minidx=minidx) 
    # tt, tv = peak.conservative_peaks(omtime, ommag, minidx=minidx, pktype="trough")

    # pspl = UnivariateSpline(pt, pv, k=5, s=0)
    # tspl = UnivariateSpline(tt, tv, k=5, s=0)
     
    # r_p = np.sqrt(pspl(int_t))
    # r_a = np.sqrt(tspl(int_t))
    # ecc = (r_p - r_a)/(r_p + r_a)
    # return np.array(ecc)


def omega_smallest(time, data, windows=range(1, 7), coord="r"):
    if coord=="r":
        oms = [omega_r_sliding(time, data, w) for w in windows]
    elif coord=="phi":
        oms = [omega_phi_sliding(time, data, w) for w in windows]
    outarr = np.zeros((len(windows), len(time)))
    for n in range(0, len(windows)):
        om = oms[n]
        omspl = UnivariateSpline(om[:, 0], om[:, 1], k=4, s=0)
        ltime = om[0, 0]
        lidx = np.argmax(time > ltime) - 1
        rtime = om[-1, 0]
        ridx = np.argmax(time > rtime)
        outarr[n, lidx:ridx] = omspl(time[lidx:ridx])
        outarr[n, 0:lidx] = 100.
        outarr[n, ridx:] = 100.
    least_om = np.zeros((len(time), 2))
    least_om[:, 0] = time
    for tidx in range(0, len(time)):
        if outarr[1, tidx] > outarr[0, tidx]:
            least_om[tidx, 1] = outarr[0, tidx]
        else:
            minom =  np.min(outarr[:, tidx])
            least_om[tidx, 1] = minom 
    goodidxs = least_om[:, 1] < 100.
    good_om = least_om[goodidxs, :]
    return good_om

def omega_phi_from_peaks(time, data, pt, pv, N, trimdelta=0.1, joined=False):
    if N<1:
        raise ValueError("Invalid number of windows N = "+str(N))
    newtime, newdat, newidx = inject_and_sort(time, data, pt, pv)
    omphN = omega_phi_do_work(newtime, newdat, newidx, 
                              pt, N, joined=joined, trimdelta=trimdelta)
    return omphN 

def omega_phi_sliding(time, data, N, trimdelta=0.1, pktype="peak",
                      refine=True, refine_pktype="peak", Npk=5,
                      polyorder=3, splorder=3, minidx=200,
                      peak_method="spline"):
    pt, pv, joined = get_peaks(time, data, pktype=pktype, refine=refine,
                               refine_pktype=refine_pktype, Npk=Npk,
                               polyorder=polyorder, splorder=splorder,
                               minidx=minidx, peak_method=peak_method)
    omphN = omega_phi_from_peaks(time, data, pt, pv, N, trimdelta=trimdelta, 
                                 joined=joined)
    return omphN

def omega_r_from_peaks(time, data, pt, N, trimdelta=0.1, joined=False):
    if N<1:
        raise ValueError("Invalid number of windows N = "+str(N))
    omrN = omega_r_do_work(time, pt, N, joined=joined, trimdelta=trimdelta)
    return omrN 
  
def omega_r_sliding(time, data, N, trimdelta=0.1, pktype="peak",
                    refine=True, refine_pktype="peak", Npk=5,
                    polyorder=3, splorder=3, minidx=200,
                    peak_method="spline"):
    """Return omega_r, calculated from peak to peak and trough to trough ranges.
       These are returned separately, as two lists, the first for the peaks
       and the second for the troughs. Each list contains a (boolean) array
       of masks which are False for extrapolated times, a list of times for which
       no interpolation has been performed, and the data fit to a spline.

       This works by computing the time-width of each peak and assigning that value to the midpoint time.

       This function can also be used to calculate the other frequencies from 
       the actual coordinate data.
    """
    if N<1:
        raise ValueError("Invalid number of windows N = "+str(N))

    pt, pv, joined = get_peaks(time, data, pktype=pktype, refine=refine,
                               refine_pktype=refine_pktype, Npk=Npk,
                               polyorder=polyorder, splorder=splorder,
                               minidx=minidx, peak_method=peak_method)
    omrN = omega_r_do_work(time, pt, N, joined=joined, trimdelta=trimdelta)
    return omrN

def omega_r_do_work(time, tpks, N, joined=False, trimdelta=0.005):
    """Compute omega_r: the reciprocal of the width of the period-to-period 
       times.
    """
    Norbits = N#%2
    if joined:
        N = 2*N
    t0s, tfs, dTs, tmids = window_ranges(tpks, N)
    omRs = 2.0*np.pi*Norbits/dTs
    tmids_trim, omRs_trim = discard_outliers(tmids, omRs, trimdelta)

    return np.column_stack((tmids, omRs))

def omega_phi_do_work(time, data, pidx, pt, N, 
                      joined=False, trimdelta=0.005):
    """Compute omega_phi: the time-average of each period.
    """
    if joined:
        N = 2*N
    t0s, tfs, dTs, tmids = window_ranges(pt, N)
    lid = pidx[:-N]
    rid = pidx[N:]
    ints = np.zeros(len(lid))
    if len(lid) != len(rid):
        raise ValueError("lid, rid had different lengths (",lid, rid,")")
    for l, r, i in zip(lid, rid, range(0, len(lid))):
        thisint = simps(data[l:(r+1)], x=time[l:(r+1)])
        ints[i] = thisint
    #ints = np.array([np.trapz(data[l:(r+1)], x=time[l:(r+1)]) for l, r in tup])
    # sums = np.array([np.sum(data[l:(r+1)]) for l, r in tup])
    # ints = []
    # for thissum, thistup in zip(sums, tup):
      # t0 = time[thistup[0]]
      # tf = time[thistup[1]+1]
      # ints.append(thissum*(tf-t0)) 
    # ints = np.array(ints)

    omphis = np.divide(ints, dTs)

    tmids_trim, omphis_trim = discard_outliers(tmids, omphis, trimdelta)
    outarr = np.column_stack((tmids, omphis))
    return outarr

def interp_ratio(numtime, numdat, denomtime, denomdat):
    num_spl = UnivariateSpline(numtime, numdat, s=0)
    return num_spl(denomtime)/denomdat

def avg_error(newt, pk, oldt, tr):
    tr_spl = UnivariateSpline(oldt, tr, s=0)
    tr_int = tr_spl(newt) 
    avg = 0.5*(tr_int+pk)
    pkdiff = pk-avg
    trdiff = tr_int-avg
    return avg, pkdiff, trdiff 


# def krph(time, data, N, Niter=1, trim=True, minidx=200):
    # """Return the omega_r/omega_phi. Data are interpolated onto the 
       # omega_r time points.
    # """
    # ph_pk, ph_tr = omega_phi_sliding(time, data, N, 
                                     # Niter=Niter,
                                     # trim=trim, minidx=minidx)
    # r_pk, r_tr = omega_r_sliding(time, data, N, Niter=Niter, trim=trim, 
                                 # minidx=minidx)

    # ph_all = omega_phi_sliding(time, data, N, Niter=Niter, trim=trim,
                               # minidx=minidx, joined=True)
    # r_all = omega_r_sliding(time, data, N, Niter=Niter, trim=trim,
                            # minidx=minidx, joined=True)
    # k_pk = r_pk[:, 1] / ph_pk[:, 1]
    # k_tr = r_tr[:, 1] / ph_tr[:, 1]
    # k_all = r_all[:, 1] / ph_all[:, 1] 

    # out = dict()
    # out["t_pk"] = ph_pk[:, 0]
    # out["t_tr"] = ph_tr[:, 0]
    # out["t_all"] = ph_all[:, 0]
    
    # out["krph_pk"] = k_pk 
    # out["krph_tr"] = k_tr 
    # out["krph_all"] = k_all

    # out["omr_pk"] = r_pk[:, 1]
    # out["omr_tr"] = r_tr[:, 1]
    # out["omr_all"] = r_all[:, 1]

    # out["omph_pk"] = ph_pk[:, 1]
    # out["omph_tr"] = ph_tr[:, 1]
    # out["omph_all"] = ph_all[:, 1]
    # return out

# def krphplot(stem, fnames, Niter=1, Nwindow=1):
    # dat=dict()
    # for name in fnames:
      # f = open(stem+name+"/Omega.pkl")
      # Omega, OmegaMag = pickle.load(f)
      # f.close()

      # out = krph(OmegaMag[:, 0], OmegaMag[:, 1], Nwindow, Niter=Niter)
      # out["ommag_t"] = OmegaMag[:, 0]
      # out["ommag"] = OmegaMag[:, 1]
      # time = hfile.get_coord("t", stem+name+"/Horizons.h5")
      # rmag = hfile.get_coord("r", stem+name+"/Horizons.h5")
      # ecc = ecc_kepler(time, rmag, out["t_all"], Niter=Niter)
      # out["ecc"] = ecc

      # dat[name] = out
    # return dat


# def ecc_kepler(time, rmag, int_t, factor=1., Niter=1, minidx=200):
    # pt, pv, pspl = peak.oscillation_peaks(time, rmag, pktype="peak",
                                          # Niter=Niter, minidx=minidx)
    # tt, tv, tspl = peak.oscillation_peaks(time, rmag, pktype="trough",
                                          # Niter=Niter, minidx=minidx)

    # r_p = pspl(int_t)
    # r_a = tspl(int_t)
    # ecc = factor*(r_p - r_a)/(r_p + r_a)
    # return np.array(ecc)


##############################################################################
#FREQUENCY SPECIFIC REDIRECTS
##############################################################################
# def integrator(time, data, dt):
    # from scipy.integrate import romb
    # from scipy.integrate import simps
    # #construct ranges to integrate over
    # ranges = [0]
    # rombpts = np.array([2.**k + 1 for k in range(0, 15)])
    # print rombpts
    
    # npts = len(data)
    # def next_range_entry(rombpts, npts, ranges):
        # shifted = rombpts - npts
        # try:
            # smallest_diff = shifted[shifted<=0][-1]
        # except IndexError:
            # return -1 
        # closest_rombpts = smallest_diff + npts
        # return closest_rombpts
    
    # nxt = 0 
    # while nxt != -1:
        # nxt = next_range_entry(rombpts, npts, ranges)
        # if nxt!=-1:
            # npts = npts - nxt
            # ranges.append(ranges[-1] + nxt)
    # ints = [] 
    # for lrmb, rrmb in zip(ranges[:-1], ranges[1:]):
        # thisint = romb(data[lrmb:rrmb], dx=dt)
        # ints.append(thisint)

    # lastl = ranges[-1]
    # if lastl != len(data):
        # lastint = simps(data[lastl:], x=time[lastl:])
        # ints.append(lastint)

    # total = np.sum(np.array(ints))
    # return total

def discard_outliers(ts, data, delta=0.1):
    """Drops frequencies which differ from their predecessors by more than
       delta.
    """
    if delta==0:
      return ts, data
    diffs = np.fabs(np.diff(data))
    #print diffs
    drops = diffs > delta
    badinds = np.where(drops)[0] + 1
    ts_trim = ts[~badinds]
    data_trim = data[~badinds]
    return ts_trim, data_trim


def trim_time_max(peaks, troughs, time):
    """Return time, trimmed so it doesn't exceed the last extremum. Useful
       to avoid using extrapolated data.
    """
    maxpeak = np.amax(peaks)
    maxtrough = np.amax(troughs)
    maxt = np.amax([maxtrough]+[maxpeak])
    critical_index = np.amin(np.where( time > maxt ))
    return time[:critical_index]

def trim_junk_time(time, data, tjunk=300):
    critical_index = np.amax(np.where( time < tjunk ))
    return time[critical_index:], data[critical_index:]

def inject_and_sort(time, data, pt, pv):
    """Insert the values pv(pt) to the numerical function data(time). 
       Return the new times, new data, and the indices of the values we 
       inserted.
    """
    inject_idxs = np.searchsorted(time, pt)
    newtime = np.insert(time, inject_idxs, pt)
    newdata = np.insert(data, inject_idxs, pv)
    newidxs = inject_idxs + np.arange(0, len(inject_idxs))
    return newtime, newdata, newidxs

def extend_first_to_second(arr1, arr2):
    """Extend (t1, d1) to include the points of (t2, d2) which are outside
       its range.
    """
    
    # arr1 = np.transpose(np.array([t1, d1]))
    # arr2 = np.transpose(np.array([t2, d2]))
    tsmall = arr1[0, 0]
    tlarge = arr1[-1, 0]
    if tsmall < arr2[0, 0] or tlarge > arr2[-1, 0]:
        raise ValueError("Inconsistent input; were arrays in wrong order?")
    
    lidx = np.argmin(np.abs(arr2[:, 0]-tsmall))
    arrleft = arr2[:lidx, :]    
    
    ridx = np.argmin(np.abs(arr2[:, 0]-tlarge))
    arrright = arr2[ridx+1:, :] #this will be empty if t1[-1]=t2[-1]
    
    lstacked = np.vstack([arrleft, arr1])
    rstacked = np.vstack([lstacked, arrright])
    return rstacked


#computes lists t0s, tfs, dTs, and tmids over windows.
def window_ranges(tpks, N):
    t0s       = np.array(tpks[:-N])
    tfs       = np.array(tpks[N:])
    dTs       = tfs-t0s 
    tmids     = t0s + dTs/2. 
    return t0s, tfs, dTs, tmids


# def extrapolate_to_N0(rng, mask, vals, deg=3):
  # nnonzero = np.sum(mask, axis=0) 
  # extrap = np.zeros(vals.shape[1]) 
  # for t_i in range(0, vals.shape[1]): 
    # thismask = mask[:, t_i]
    # thisvals = vals[:, t_i]
    # thisdeg = deg
    # extrap[t_i] = thisvals[0]
    # #Fit becomes badly conditioned if very few windows are unmasked
    # #while (2*thisdeg >= nnonzero[t_i]) and (thisdeg>0):
    # #  thisdeg = thisdeg - 1
    # #if 0==thisdeg:
    # #  extrap[t_i] = thisvals[0]
    # #else:
    # #  z = np.polyfit(rng, thisvals, w=thismask, deg=thisdeg) 
    # #  p = np.poly1d(z) 
    # #  extrap[t_i]=p(0)
    # #print "tstep = ", t_i, "nnonzero = ", nnonzero[t_i], "thisdeg = ", thisdeg, "extrap = ", extrap[t_i]
  # return extrap

# def make_om_data(omtype, time, data, rng):
  # if omtype=="r":
    # func = omega_r_sliding 
  # elif omtype=="phi":
    # func = omega_ph_sliding 
  # else:
    # raise NotImplementedError("'r' and 'phi' only supported options")
  # omarr = np.zeros([len(rng), len(time), 2])
  # maskarr = np.zeros([len(rng), len(time), 2], dtype=bool)
  # for Ni in range(0, len(rng)):
    # peakmask, peakdata, troughmask, troughdata = func(time, data, rng[Ni])
    # maskarr[Ni, :, 0] = np.copy(peakmask)
    # omarr[Ni, :, 0] = np.copy(peakdata)
    # maskarr[Ni, :, 1] = np.copy(troughmask)
    # omarr[Ni, :, 1] = np.copy(troughdata)
  # return maskarr, omarr

# #Return data vs. N at t=t0
# def extrapolation_plot(omtype, t0, time, data, rng=range(1, 30, 1)):
  # maskarr, omarr = make_om_data(omtype, time, data, rng)
  # return [maskarr[:, t0, :], omarr[:, t0, :]] 

# #The main function to be called from outside
# def extrapolated_omega(omtype, time, data, rng=range(1, 30, 1)):
  # maskarr, omarr = make_om_data(omtype, time, data, rng)
  # pextrap = extrapolate_to_N0(rng, maskarr[:, :, 0], omarr[:, :, 0], deg=3)
  # textrap = extrapolate_to_N0(rng, maskarr[:, :, 1], omarr[:, :, 1], deg=3)
  # both = np.column_stack([pextrap, textrap])
  # extrap = (pextrap+textrap)/2.
  # terr   = np.fabs(textrap-extrap)
  # perr   = np.fabs(pextrap-extrap)
  # err    = (terr+perr)/2.
  # return extrap, err


