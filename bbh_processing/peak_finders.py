"""This collects functions to find peaks of orbital time series.
"""
from scipy.interpolate import UnivariateSpline, lagrange
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.signal import savgol_coeffs

def drop_duplicate_columns(df):
    """Deletes duplicate columns from df.
    """
    return df.T.groupby(level=0).first().T

def apply_func(row, newname, func, col1):
    result = func(row[col1])
    return pd.Series({newname: result})

def add_omegamag(df, func, col1):
    assert df is not None, "Passed in empty df!"
    def to_apply(row):
        time, omegamag = func(row[col1])
        return pd.Series({"omtime": time, "omega_total": omegamag})
    listdf = df.apply(to_apply, axis=1)
    concat = pd.concat([df, listdf], axis=1)
    return drop_duplicate_columns(concat)



def apply_func2(row, newname, func, col1, col2):
    result = func(row[col1], row[col2])
    return pd.Series({newname: result})

def add_column(df, newname, func, intag):
    assert df is not None, "Passed in empty df!"

    if len(intag)==1:
      def to_apply(row):
          return apply_func(row, newname, func, intag[0])
    if len(intag)==2:
      def to_apply(row):
          return apply_func2(row, newname, func, intag[0], intag[1])

    listdf = df.apply(to_apply, axis=1)
    concat = pd.concat([df, listdf], axis=1)
    return drop_duplicate_columns(concat)

def spline_quantities(df, col1="Time", col2="r_mag"):
    """This is the thing to call from outside.
    """
    assert df is not None, "Passed in empty df!"
    def apply_peak(row):
        return spline_quantities_series(row[col1], row[col2], k=4, s=0)
    listdf = df.apply(apply_peak, axis=1)
    concat = pd.concat([df, listdf], axis=1)
    return drop_duplicate_columns(concat)

def add_tags(df):
    """Adds identifying tags to the runs.
    """

    def make_string(row):
      eccstr = '%.2f' % row.EccN
      incstr = '%.f' % row.IncA
      chistr = '%.2f' % row.ChiA
      dstr = '%.f' % row.D0
      return pd.Series({"tag":"EccN_"+eccstr+"IncA_"+incstr+"ChiA_"+chistr+"D_"+dstr})
    tagdf = df.apply(make_string, axis=1)
    concat = pd.concat([df, tagdf], axis=1)
    return drop_duplicate_columns(concat)

def trim_junk(to_trim, junkthresh=100, keepfirst=True):
    #mask = np.ones(len(to_trim), dtype=bool)
    mask = to_trim>junkthresh
    if keepfirst==False:
        return_array = to_trim[mask]
    else:
        temp = to_trim[mask]
        return_array = np.zeros(len(temp)+1)
        return_array[0] = 0.
        return_array[1:] = temp[:]
    return return_array
    #diff = np.diff(to_trim)

    #junk_regime = True
    #i = 0
    #maxlen = len(to_trim)
    #while junk_regime and i < maxlen:
    #    if diff[i] <= junkthresh:
    #        mask[i] = False
    #    else:
    #        junk_regime = False
    #    i = i + 1



###############################################################################
#PEAK FINDERS
# Functions that do the brute numerical work of finding peaks.
###############################################################################

def joined_peaks(time, data, refine=True, refine_pktype="peak", Npk=5, 
                 polyorder=3, splorder=3, minidx=200, peak_method="spline",
                 Nref=8):
    """Return the peaks and troughs, joined together in order.
    """
    ets = []
    evs = []
    for pktype in ["peak", "trough"]:
        et, ev = conservative_peaks(time, data, pktype=pktype, refine=refine,
                                    refine_pktype=refine_pktype,
                                    Npk=Npk, polyorder=polyorder, 
                                    splorder=splorder, minidx=minidx, 
                                    peak_method=peak_method, Nref=Nref)
        ets.append(et)
        evs.append(ev)
    #Trim to the length of the shortest
    allt = np.concatenate(ets)
    allv = np.concatenate(evs)
    joined = np.transpose(np.array([allt, allv]))
    sort_arr = joined[joined[:, 0].argsort()]
    return sort_arr[:, 0], sort_arr[:, 1]
   
def conservative_peaks(time, indata, pktype="peak", refine=True,
                       refine_pktype="peak", Npk=5, polyorder=3, splorder=3, 
                       minidx=200, peak_method="spline", Nref=8,
                       returnrefined=False):
    """Given a time-series (time, data) which samples a function F(t) with an 
       oscillatory ('conservative', though possibly chirping) O(t) plus a 
       monotonic ('dissipative') component M(t) (so F(t) = O(t) + M(t)), 
       estimate the class 'pktype' extrema of O(t).

       This function first computes the class-'refine_pktype' extrema of O(t).  
       Note this class may be different from those of 'pktype' since we 
       generally find that periastra produce better estimates of the underlying
       trends than do apoastra. Generally this function is applied to the 
       mean frequency of a BBH inspiral, in which case peaks correspond to
       periastron passages. If applied to for example radial data the reverse
       is true and refine_pktype should generally be set to "trough".
       
       The options Npk, polyorder, and peak_method
       control the behaviour of the peak finder, which in general works by first
       finding guess points in the data which are not monotonic in their 
       surroundings. Guess peaks which are "too close together" are smoothed
       over using spline interpolation.

       The guess points are then improved upon by extremizing the derivative
       of some interpolant. The exact nature of the interpolant depends on the
       choice of 'peak_method', but by default is an order-4 spline
       fitted to ranges of Npk points to the left and right of the data 
       (this appears to produce best results for curves from
       BBH inspirals which typically exhibit significant asymmetry about
       their extrema).

       We then extremize the derivative of that spline and return the root
       closest to the guess point. Note that if we are using spline 
       interpolation the option 'polyorder' has no effect. For the other 
       methods - notably 'polyfit', which interpolates using a single polynomial
       rather than a piecewise spline - 'polyorder' sets the order of the 
       interpolant.

       Next it computes 'omega_phi', which is an estimate of the time-weighted
       rolling mean of the overall function (and thus of M(t)). Specifically, 
       we take M_t = Int(F(t))/delta(t), where Int(F(t)) is the integral of 
       F(t) between the extrema of F(t) and delta(t) is the time interval 
       between them. The time points of M_t are taken to be halfway between
       the peaks.

       We now approximate M(t) as the order-'splorder' spline interpolating
       M_t to the original time-series of F(t). This will normally involve
       some extrapolation at the boundaries, since the first (last) value of
       omega_phi will occur somewhat to the right (left) of F(t).

       With M(t) computed we now take O(t) = F(t) - M(t) and repeat
       the peak finding procedure as defined above. We then add M(t) to the
       values (not the times) of these new peaks in order to produce our
       final result.

       Minidx sets a minimum index prior to which peaks are ignored. This is 
       useful for excluding spurious peaks from junk radiation. The default
       value of 200 avoids the junk well in practice.

       For convenience, this function can also be instructed to return the 
       estimated peaks of F(t) without refinement, by setting the option
       'refine' to False.
    """

    #Compute the peaks of the original function. If we are not refining, 
    #return them.
    data = np.copy(indata)
    
    #omph_spl = UnivariateSpline(omph[:, 0], omph[:, 1], k=splorder, s=0)
    for i in range(0, Nref): 
        if not refine:
          refine_pktype = pktype
        pt, pv, pspl = envelope(time, data, envtype=refine_pktype, Npk=Npk,
                                polyorder=polyorder, splorder=splorder,
                                minidx=minidx, method=peak_method) 
        if not refine:
            return [pt, pv]
        data = data - pspl(time)
        # plt.plot(time, data)
        # plt.ylim((-0.01, 0.01))
        # plt.xlim((23000, 23300))
        # plt.show()
        #cons_spl = UnivariateSpline(time, cons_data, k=4, s=0)

        #Subtract the time-weighted mean (omega_phi).
        # omph = eqfreq.omega_phi_from_peaks(time, data, pt, pv, 1)  
        # omph_interp = omph_spl(time)
        #conservative_data = data - omph_interp
    #Now compute the peaks a second time. These estimate the peaks of the
    #underlying conservative trend in the data.
    ct, cv, cspl = envelope(time, data, envtype=pktype, Npk=Npk,
                            polyorder=polyorder, splorder=splorder, 
                            minidx=minidx, method=peak_method)
    spl0 = UnivariateSpline(time, indata, k=4, s=0)
    cv_refined = spl0(ct)
    returnlist = [ct, cv_refined]
    if returnrefined:
        returnlist.append(data)
        returnlist.append(cv)
    return returnlist



def oscillation_peaks(time, data, pktype="peak", Niter=1,
                      Npk=5, polyorder=2, splorder=3, minidx=200, 
                      return_spline=True, maintain_length=True):
    """Give a time-series (time, data) which is assumed to be sampled from a 
       monotonic plus an oscillating function, attempt to find the peaks of
       the oscillating part only. These will generally be shifted somewhat
       from the extrema of the total.

       To do this, we first find the peaks (troughs) of the total using the 
       polynomial method and construct an envelope function from them. This
       is the zeroth iteration q0.
       
       Next, we interpret q0 as an estimate of the monotonic part of the 
       function. Thus the desired peaks (troughs) are extrema of [data-q0], not 
       of data alone. We thus construct this difference, extremize it once 
       again, and recompute the envelope, yielding the first iteration q1.

       This process is then repeated Niter times (default Niter=1, returning
       q1).

       Sometimes, the refined enveloped has peaks which were wholly absent
       from the original. If maintain_length is True, these are excluded from
       the output.
    """
    def do_return(qit, qiv, qispl):
        if return_spline:
          return [qit, qiv, qispl]
        else:
          return [qit, qiv]


    def refine(oldtime, olddat, q0spl, 
               Npk=Npk, polyorder=polyorder, splorder=splorder,
               minidx=minidx, old_length=-1): 
        qinterp = q0spl(oldtime)
        shift = olddat - qinterp
        
        #here we assume einterp is monotonic
        q1t, sval, shift_spl = envelope(oldtime, shift, envtype=pktype,
                                        Npk=Npk, polyorder=polyorder,
                                        splorder=splorder, minidx=minidx)
        
        if old_length != -1:
          if len(q1t) > old_length:
            q1t = q1t[:old_length]
            sval = sval[:old_length]

        q1val = sval + q0spl(q1t)
        q1spl = UnivariateSpline(q1t, q1val, s=0, k=splorder)
        return [q1t, q1val, q1spl]

    qit, qiv, qispl = envelope(time, data, envtype=pktype, 
                               Npk=Npk, polyorder=polyorder,
                              splorder=splorder, minidx=minidx)
    if 0==Niter:
        return do_return(qit, qiv, qispl)
    
    trunc_time = time[minidx:]
    trunc_dat = data[minidx:]
    for i in range(0, Niter):
        if not maintain_length:
          oldlength = -1
        else:
          oldlength = len(qit)
        qit, qival, qispl = refine(trunc_time, trunc_dat, qispl, 
                                   Npk=Npk, polyorder=polyorder, 
                                   splorder=splorder, minidx=minidx,
                                   old_length=oldlength)
    
    return do_return(qit, qival, qispl)

def envelope(time, data, envtype="peak",
             Npk=10, polyorder=2, splorder=3, minidx=200, 
             method="spline"):
    """Given a time-series (time, data), finds all the extrema of the class
       specified by 'envtype'. Extrema are computed by finding the (assumed 
       unique) root of an order 'polyorder' polynomial fit to Npk points to the
       left and right of a guess value found by comparing adjacent points in 
       the time series. The first 'minidx' points are ignored (this allows
       junk-contaminated data to be excluded).
       The returned list contains, in order:
         1. The time-values of the peaks.
         2. The data-values of the peaks.
         3. An order-'splorder' UnivariateSpline passing through all the 
            peaks.
    """
  
    pidx, pvg, tidx, tvg = turning_points(data, minidx=minidx)
    #ptg = time[pidx]
    #ttg = time[tidx]
    #ptg, pvg, ttg, tvg = peaks_from_splines(time, data, minidx=minidx)
    if envtype=="peak":
        #etg, evg = ptg, pvg
        eidx, ev = pidx, pvg
    elif envtype=="trough":
        #etg, evg = ttg, tvg
        eidx, ev = tidx, tvg
    else:
        raise ValueError("Invalid envtype " + envtype)
    etimes, evals = peaks_from_polyfit_work(time, data, eidx, Npk, polyorder,
                                            method=method, input_type="index")
    env_spline = UnivariateSpline(etimes, evals, s=0, k=splorder)
    # plt.plot(time, data)
    # plt.plot(etimes, env_spline(etimes))
    # plt.show()
    return [etimes, evals, env_spline]


# def peaks_from_polyfit(time, data, N=5, order=2, minidx=200):
    # """First finds approximate locations of peaks using exact data. Next applies
       # a polynomial interpolant in an N-point range to either side of the
       # peak.
    # """
    # pidx, pv, tidx, tv = turning_points(data, minidx=minidx)
    # ptimes, pvals = peaks_from_polyfit_work(time, data, pidx, N, order)
    # ttimes, tvals = peaks_from_polyfit_work(time, data, tidx, N, order)
    # return ptimes, pvals, ttimes, tvals

def peaks_from_polyfit_work(time, data, pt, N, order, 
                            input_type="time", method="spline"):
    """This does the heavy lifting for 'peaks_from_polyfit'. 
    """
    npts = len(pt) 
    extrema_times = np.zeros(npts)
    extrema_vals = np.zeros(npts)
    i = 0
    for t in pt:
        if input_type=="time":
            idx = np.argmin(np.abs(time-t))
        elif input_type=="index":
            idx = t
        else:
            raise ValueError("Invalid input_type")
        #lidx = idx-N
        ridx = idx+N
        if ridx >= len(time):
            ridx = len(time)-1
        tdiff = time[ridx] - time[idx]
        lidx = np.argmin(np.abs(time - (time[idx] - tdiff)))
       
        #print lidx, ridx
        if lidx <= 0: #or ridx >= len(time):
            extrema_times[i] = time[idx]
            extrema_vals[i] = data[idx]
        else:
            fit_times = time[lidx:ridx]
            fit_vals = data[lidx:ridx]
            et, ev = fit_extremum(fit_times, fit_vals, method, 
                                  time[idx],
                                  order=order, spacing=0.5)
            extrema_times[i] = et
            extrema_vals[i] = ev
            i = i + 1
    return extrema_times, extrema_vals

def fit_extremum(times, vals, method, closest, order, spacing=0.5, window=7):
    if method=="polyfit" or method=="lagrange":
        if method=="polyfit":
            interp = np.poly1d(np.polyfit(times, vals, order))
        elif method=="lagrange": #warning - this is unstable
            interp = lagrange(times, vals)  
        deriv = np.polyder(interp)
        roots = np.roots(deriv)
        extremum_idx = np.argmin(np.abs(roots - closest))
        et = roots[extremum_idx]
        ev = interp(et)
        return et, ev

    elif method=="spline":
        spl = UnivariateSpline(times, vals, k=4, s=0)
        extrema = spl.derivative().roots()
        try:
            extremum_idx = np.argmin(np.abs(extrema - closest))
        except ValueError:
            raise ValueError("Found no extrema in range : ", times)
        et = extrema[extremum_idx]
        ev = spl(et)
        return et, ev



    elif method=="SG":
        #fit = savgol_filter(vals, window, order, deriv=0, delta=spacing)
        deriv = savgol_filter(vals, window, order, deriv=1, delta=spacing)
        goodtimes = times[window:-(window//2)]
        goodvals = deriv[window:-(window//2)]
        spl = UnivariateSpline(goodtimes, goodvals, k=3, s=0)
        extrema = spl.roots()
        extremum_idx = np.argmin(np.abs(extrema - closest))
        et = extrema[extremum_idx]
        ev = spl(et)
        return et, ev
    else:
        raise NotImplementedError("Method " + method + " not implemented.")


# def peaks_from_polyfit_work(time, data, pt, N, order):
    # """This does the heavy lifting for 'peaks_from_polyfit'. 
    # """
    # npts = len(pt) 
    # extrema_times = np.zeros(npts)
    # extrema_vals = np.zeros(npts)
    # i = 0
    # for t in pt:
        # idx = np.argmin(np.abs(time-t))
        # lidx = idx-N
        # ridx = idx+N
        # if lidx <= 0 or ridx >= len(time):
            # extrema_times[i] = time[idx]
            # extrema_vals[i] = data[idx]
        # else:
            # fit_times = time[lidx:ridx]
            # fit_vals = data[lidx:ridx]
            # interp = np.poly1d(np.polyfit(fit_times, fit_vals, order))
            # extrema = np.roots(np.polyder(interp))
            # extremum_idx = np.argmin(np.abs(extrema - time[idx]))
            # extrema_times[i] = extrema[extremum_idx]
            # extrema_vals[i] = interp(extrema_times[i])
            # i = i + 1
    # return extrema_times, extrema_vals

def peaks_from_splines(time, data, k=4, s=0, minidx=200):
    """Applies an order-k spline to the data and finds the roots of its
       derivative. These are the times at which the function
       peaks. Also returns the value of the spline at these times.
       Optionally, returns the spline itself.
    """
    if k < 4:
        print "Warning: k must be at least 4; using 4 instead of", k
        k = 4
    spl = UnivariateSpline(time[minidx:], data[minidx:], k=k, s=s)
    peaktimes = spl.derivative().roots()
    #peaktimes = trim_junk(peaktimes, junkthresh=100, keepfirst=True)
    peakvals = spl(peaktimes)
    return_list = peaks_and_troughs(peaktimes, peakvals)
    #return_list.append(spl)
    return return_list


def interpolate_peak(peaktime, peakval, time, maxk=3):
    """Interpolate a peak to time 'time'. Does different things for
       different amounts of data.
    """
    plen = len(peaktime)
    assert plen == len(peakval), "Inconsistent lengths."
    if plen == 0:
        raise ValueError("No data in peaktime!")
    if plen <= 1:
        interpval = np.zeros(len(time))
        interpval.fill(peakval[0])
    else:
        myk = plen - 1
        if myk > maxk:
            myk = maxk
        spl = interpolate.UnivariateSpline(peaktime, peakval, k=myk, s=0)
        interpval = spl(time)
    if type(interpval)==np.float64:
        interpval = np.array([interpval])
    return interpval

def compute_orbital_quantities(peaktimes, peakvals, troughtimes, troughvals,
                               norbits, time, data):
    """Compute r_a, r_p, eccentricity, semi-major axis, semi-latus rectum.
    """
    if 0 == norbits:
        return None, None, None, None, None
    else:
        r_a = interpolate_peak(peaktimes, peakvals, time, maxk=3)
        r_p = interpolate_peak(troughtimes, troughvals, time, maxk=3)
        replace_inds = np.where(r_p >= r_a)[0]
        r_p[replace_inds] = np.copy(r_a[replace_inds])
        ecc = (r_a - r_p) / (r_p + r_a)
        semi_maj = (r_a + r_p)/2.
        rho = 2. / (1./r_a + 1./r_p)
    return r_a, r_p, ecc, semi_maj, rho

def join_plunge(series, peaktimes, time, data):
    """Detect the plunge; use the actual r data there rather than peak
    extrapolations.
    """
    deriv = np.gradient(np.gradient(data, np.gradient(time)),np.gradient(time))
    plungeinds = np.where(deriv<-0.003)[0]
    if len(plungeinds) == 0:
        return series
    firstplungeind = np.amin(plungeinds)

    series[firstplungeind:] = data[firstplungeind:]
    return series



def spline_quantities_series(time, data, k=4, s=0):
    """Constructs a pandas series of various useful quantities.
    """
    peaktimes, peakvals, troughtimes, troughvals, spl = peaks_from_splines(
                                                              time, data, k, s)
    norbits = len(peaktimes) - 1
    if len(troughtimes) > len(peaktimes):
        norbits += 0.5
    r_a, r_p, ecc, semi_maj, rho = compute_orbital_quantities(
          peaktimes, peakvals, troughtimes, troughvals, norbits, time, data)

    return pd.Series(
                      {'peaktimes': peaktimes,
                       'peakvals': peakvals,
                       'troughtimes': troughtimes,
                       'troughvals': troughvals,
                       'r_spl': spl,
                       'r_p': r_p,
                       'r_a': r_a,
                       'ecc': ecc,
                       'semi_maj': semi_maj,
                       'rho': rho,
                       'norbits': norbits
                      }
                    )


def simple_peak_test(data, pnt):
    """Returns true if data[pt] is a peak (i.e. its neighbours are both either
    greater or lesser than it).
    """
    #maximum
    if data[pnt] > data[pnt-1]:
        if data[pnt] > data[pnt+1]:
            return True
    #minimum
    if data[pnt] < data[pnt-1]:
        if data[pnt] < data[pnt+1]:
            return True
    #neither
    return False

def find_bumps(ev, thresh=0.0001):
    """Returns ranges spanning bumps in the data whose height is beneath a 
       certain threshold. The convoluted logic allows the range to expand
       to include possibly-multiple too-close points.
    """
    rgs = []
    #for i in range(0, len(ev)-1):
    i = 0
    while i < len(ev)-1:
        rrng = i 
        if np.abs(ev[i] - ev[i+1]) < thresh:
            lrng = i
            off = 2
            gate = True
            while gate:
                if i+off >= len(ev):
                    gate=False
                elif np.abs(ev[i+off-1] - ev[i+off]) > thresh:
                    gate=False
                else:
                    off = off + 1
            rrng = i+off-1
            rgs.append((lrng, rrng))
        i = rrng + 1

            # if np.abs(ev[i+1] - ev[i+2]) < thresh:
                # rgs.append((i, i+2))
    return rgs

def smooth_bumps(data, rgs, turns, d=5, ls=2, rs=3):
    idx = np.linspace(0, len(data), len(data), endpoint=False, dtype=np.int)
    newdata = np.copy(data)
    for rg in rgs:
        il = turns[rg[0]] - ls
        ir = turns[rg[1]] + rs
        if il > 0:
          if ir < len(newdata):
            lrng = idx[0:il]
            rrng = idx[ir:]
            good = np.concatenate((lrng, rrng))
            spl = UnivariateSpline(good, data[good], s=0, k=3)
            mrng = range(il, ir) #points to interpolate
            newdata[mrng] = spl(idx[mrng]) 
    return newdata


def get_turns(data, minidx=200):
    myrange = range(minidx, len(data)-1)
    turns = [pt for pt in myrange if simple_peak_test(data, pt)]
    turn_vals = [data[pt] for pt in turns]
    return turns, turn_vals

def turning_points(data, joined=False, minidx=200, thresh=0.001,
                   return_smooth=False):
    """Simple way to find turning points.
    """
    turns, turn_vals = get_turns(data, minidx=minidx)
    rgs = find_bumps(turn_vals, thresh=thresh) 
    newdata = smooth_bumps(data, rgs, turns)
    turns, turn_vals = get_turns(newdata, minidx=minidx)
    if joined:
        peaklist = [turns, turn_vals]
    else:
        peaklist = peaks_and_troughs(turns, turn_vals)
    if return_smooth:
        return peaklist, newdata
    return peaklist

def peaks_and_troughs(times, vals):
    """Sort data into peaks and troughs.
    """
    peaktimes = times[::2]
    troughtimes = times[1::2]
    troughvals = vals[1::2]
    peakvals = vals[::2]
    if sum(troughvals) > sum(peakvals):
        peakvals, troughvals = troughvals, peakvals
        peaktimes, troughtimes = troughtimes, peaktimes

    return [peaktimes, peakvals, troughtimes, troughvals]

#This import needs to be at the end of the file to avoid circular dependencies
#(many of the functions in equatorial_freq_methods use the peak finders)
import bbh_processing.equatorial_freq_methods as eqfreq 

