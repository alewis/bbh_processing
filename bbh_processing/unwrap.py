import numpy as np
import peak_finders as peaks
import matplotlib.pyplot as plt
from equatorial_freq_methods import window_ranges
from scipy.interpolate import UnivariateSpline
from itertools import izip as izip


def omega_theta(thtime, theta, pktime, pkdata, N=1, minidx=200, refine=True):
    #Njoined = N*2
    #Njoined = N
    pkt, pkv = peaks.conservative_peaks(pktime, pkdata, refine=True, 
                                        minidx=minidx)
    
    chi = monotonic_acos(thtime, theta, tnew=pkt, unwrap=True)
    t0s, tfs, dT, tmids = window_ranges(pkt, N)

    dChi = np.array(chi[N:] - chi[:-N])
    omths = dChi / dT
    outarr = np.column_stack((tmids, omths))
    return outarr

def omega_theta_from_peaks(thtime, theta, pkt, N=1):
    t0s, tfs, dT, tmids = window_ranges(pkt, N)

    chi = monotonic_acos(thtime, theta, tnew=pkt, unwrap=True)
    dChi = np.array(chi[N:] - chi[:-N])
    omths = dChi / dT
    outarr = np.column_stack((tmids, omths))
    return outarr

def piecewise_acos(data, ascending):
    chi = np.arccos(data)
    if ascending:
        chi = 2*np.pi - chi
    return chi

def ampoffset(time, data):
    pt, pv = peaks.conservative_peaks(time, data, pktype="peak")
    tt, tv = peaks.conservative_peaks(time, data, pktype="trough")
    tspl = UnivariateSpline(tt, tv, k=4, s=0)
    pspl = UnivariateSpline(pt, pv, k=4, s=0)
    off = tspl(time) + 0.5*( pspl(time) - tspl(time) )
    amp = 0.5*( pspl(time) - tspl(time) )
    return amp, off 

def truefalse(initial, length):
    array = np.full(length, initial, dtype=np.bool)
    array[1::2] = not initial
    return array

def monotonic_acos(time, data, tnew=None, unwrap=True, minidx=200, Niter=0):
    """Identifies increasing and decreasing branches of 'data' by finding 
       extrema
       and comparing their heights (e.g. the range between a minimum [maximum]
       and a maximum [minimum] is increasing [decreasing]). Then, applies
       the following piecewise function:
          y = acos(data) {decreasing branch}
          y = 2Pi - acos(data) {increasing branch}
       This produces a monotonically increasing phase up to jump discontinuities
       at the increasing-decreasing transition points (i.e. at the branch cuts
       of acos). If unwrap is True, this function will add factors of 2Pi to 
       the jump discontinuities in order to produce a smooth curve.

       If tnew is not None, the above function is evaluated at the points
       data(tnew) only, which are interpolated from the original time series
       'time' using a third order spline. The arguments 'time' and 'tnew'
       must both be supplied if either one is.
       
       Typically tnew will be the radial extremum times found from omega.
    """
    pidx, pv, tidx, tv = peaks.turning_points(data, minidx=minidx) 
    amp, off = ampoffset(time, data)
    dscale = (data - off) / amp
    et, ev = peaks.joined_peaks(time, data)
    
    nrngs = len(et)+1
    npi = np.zeros(nrngs)
    if tidx[0] < pidx[0]:   #initially descending 
        #ascending = [0, 1, 0, 1, 0, 1...]
        #npi = [0, 0, 1, 1, 2, 2, 3, 3...]
        ascending = truefalse(False, nrngs)
        npi[::2] = np.arange(0, len(npi[::2]))
        npi[1::2] = np.arange(0, len(npi[1::2]))
    elif pidx[0] < tidx[0]: #initially ascending
        #ascending = [1, 0, 1, 0, 1...]
        #npi = [0, 1, 1, 2, 2, 3, 3...]
        ascending = truefalse(True, nrngs)
        npi[1::2] = np.arange(1, len(npi[1::2]) + 1)
        npi[2::2] = np.arange(1, len(npi[2::2]) + 1)
    else:
        raise ValueError("First extremum identified as both trough and peak.")
    pzip = zip(et, ascending, npi)
    if tnew is not None:
        spl = UnivariateSpline(time, dscale, s=0, k=4)
        dnew = spl(tnew)
    else:
        dnew = dscale
        tnew = time
    
    chi = np.zeros(len(tnew))
    lidx = 0
    breaknow = False
    for t_r, ascend, thisnpi in pzip:
        if t_r > tnew[-1]:
            ridx = len(tnew)
            breaknow = True
        else:
            ridx = np.argmin(tnew <= t_r)
        chi[lidx:ridx] = piecewise_acos(dnew[lidx:ridx], ascend)
        if unwrap:
            chi[lidx:ridx] += 2.0*np.pi*thisnpi
        if breaknow:
            break
        lidx = ridx
    return chi 
        
        
def unwrap_acos(times, data, peak_times): 
    unwrapped = np.copy(data)
    bin_idx = np.digitize(peak_times, times) 
    jump_idx, n_pi = np.unique(bin_idx, return_counts=True)
    phase_jumps = n_pi * 2.0 * np.pi
    
    for jidx, jump in izip(jump_idx, phase_jumps):
        unwrapped[jidx:] += jump
    return unwrapped


def flip(data, junk_inds=0, shift=True):
    """Reflect descending parts of an oscillating function such that they increase.
    """
    pk_idx, pk_v, tr_idx, tr_v = peaks.turning_points(data, minidx=junk_inds)
    for a in [pk_idx, pk_v, tr_idx, tr_v]:
        if len(a) <= 0:
            return data
        #assert len(a) > 0, "No peaks found!"
    #we always want to flip between peaks and troughs
    initially_descending = True
   
    if pk_idx[0] == 0:
        initially_descending = True
    elif tr_idx[0] == 0:
        initially_descending = False
    elif pk_idx[0] < tr_idx[0]:
        initially_descending = False
    elif pk_idx[0] > tr_idx[0]: 
        initially_descending = True 
        pk_idx = [0] + pk_idx
        pk_v = [data[0]] + pk_v
    else:
        raise ValueError("Identified same point as peak and trough")
   
    last = len(data)
    if pk_idx[-1] != last and tr_idx[-1] != last:
        if pk_idx[-1] > tr_idx[-1]: #last extremum is a peak so we end descending
            tr_idx = tr_idx + [last]
    if tr_idx[0] == 0:
        tr_idx = tr_idx[1:]
    ranges = zip(pk_idx, tr_idx, pk_v)
    output = np.copy(data) 
    for l, r, this_max in ranges:
          output[l:r] = 2*this_max - data[l:r]
        
    ##account for the initial phase 
    if initially_descending and shift:
        output -= 2*output[0]
    return output


def unwrap(data, thresh=0.2, initial_fraction=1., scale=1., junk_inds=200):
    """Produce a continuous 'unwrapped phase' which increases monotonically and accretes
       by 2Pi with every period in the data.
       Step 1: make monotonic by reflecting around the y=ymax line
               at every point where the derivative changes sign.
            2: This produces a set of discontinous curves, one per cycle. Each
               should represent an increase by scale. Thus, multiply every range
               by scale*(max - min).
            3: Make the curve continuous by
               shifting each cycle up such that it joins with its predecessor.
       Params: x'thresh' is a minimum discontinuity required to identify a cycle.
                 Useful for ignoring junk.
               x'starts_ascending' should be set to True if the data is initially
                 ascending (in principle this could be detected automatically).
                 If your results seem to be correct up to a y-reflection for each
                 period, try flipping this value.
               x'intial_fraction' is the fraction of the relevant phase the orbit
                 begins at. For r-data, for example, we typically begin at 
                 apastron and this should be set to 0.5. If your initial cycle
                 seems distorted compared to the others in a systematic (i.e.
                 non-junk like way) you probably need to adjust this.
    """
    flipped = flip(data, junk_inds=junk_inds)
    pk_idx, pk_v, tr_idx, tr_v = peaks.turning_points(flipped, minidx=junk_inds)
    for a in [pk_idx, pk_v, tr_idx, tr_v]:
        if len(a) <= 0:
            return flipped, flipped
    output = np.copy(flipped)
    if pk_idx[0] == 0:
        zipvals = zip(pk_idx[1:], pk_v[1:])
    else:
        zipvals = zip(pk_idx, pk_v)
    for idx, val in zipvals:
        output[idx+1:] += val #+ central 
    return output, flipped
    
    # diff = np.abs(np.diff(flipped))
    # idx = np.arange(len(flipped))
    # flip_idx_raw=idx[diff>thresh]+1 
    # flip_idx = np.zeros(len(flip_idx_raw)+1)
    # flip_idx[1:] = flip_idx_raw

    # summands = np.ones(len(flip_idx))
    # summands[0] = 0
    # summands *= scale
    # summands[1]*= initial_fraction
    # summands = np.cumsum(summands)
    # print summands
    # output = np.copy(flipped)
   
    # frac = initial_fraction
    # for i in range(0, len(flip_idx)-1):
        # l = flip_idx[i]
        # r = flip_idx[i+1]
        # these_data = flipped[l:r]
        # shifted = these_data - these_data[0]

        # this_scale = 2*np.pi*frac/(np.amax(shifted))
        # frac = 1.
        # output[l:r] = shifted * this_scale
        # output[l:r] += summands[i]
    # #final cycle
    # last = flip_idx[-1]
    # last_data = flipped[last:]
    # shift = last_data - np.amin(last_data) 
    # output[last:] = shift*this_scale
    # output[last:] += summands[-1]
