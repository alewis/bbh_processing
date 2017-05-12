""" Tools to create a 'Mino time' phase enabling accurate integration over
    radial periods.
"""
import numpy as np
import peak_finders as peaks
from scipy.optimize import leastsq
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import warnings
from unwrap import unwrap
from cosine_power_series import CosinePowerSeries 
import bbh_processing.equatorial_freq_methods as freq 

def power_series(time, data, order=6, divvy_data=None, coef=0.6): 
    """ Compute a list of power series fits and times they are fit to of
        order 'order'. The time ranges are found by computing the minima
        of 'divvy_data' and taking ranges of coef*(time_between_minima) around
        each one. divvy_coef should generally be greater than 0.5 so that one
        need not evaluate fits exactly at the edges, but not too large to keep
        the fits accurate at low orders.
    """
    pieces = divvy_coef(time, divvy_data, coef=coef)
    fits = [CosinePowerSeries(time[l:r], data[l:r], order) for l, r in pieces]
    times = [time[l:r] for l, r in pieces]
    inner_pieces = divvy_coef(time, divvy_data, coef=0.5)
    return [times, fits, pieces, inner_pieces] 

def omega_theta_chi(time, theta, peak_data, fit_order=6, peak_coef=0.6,
                    fit_func=power_series):
    """ The "main function" to compute omega_theta. Finds the peaks of 'peak
        data' and takes ranges of coef*(time_between_minima) around each one.
        The function 'fit_funct' is then called, which applies a fitting 
        function to each range. A phase is computed from the resulting fit 
        function, which is joined together and interpolated. Finally, 
        omega_theta is the ratio of the phase and time differences between
        successive peak_data minima.
    """
    times, fits, pieces, inner = power_series(time, theta, order=fit_order, 
                                              divvy_data=peak_data, 
                                              coef=peak_coef)
    phases = offset_phases(times, fits, time, peak_data)
    joined = joined_phase(times, phases, time, peak_data)
    spl = interp.UnivariateSpline(time, joined, k=5, s=0) 
    #pt, pv, tt, tv = peaks.peaks_from_polyfit(time, peak_data)
    et, ev = peaks.joined_peaks(time, peak_data)
    om_th=[]
    this_t=[]
    for lt, rt in zip(et[:-2], et[2:]):
        deltaChi = spl(rt) - spl(lt)
        deltaT = rt - lt
        om_th.append(deltaChi / deltaT)
        this_t.append(lt + 0.5*deltaT)
    return (this_t, om_th, joined)

def phases(times, fits):
    """ Calculate the phases of the power series computed above. 
        This does not join the phases together (i.e. there will be 
        discontinuities between each fit.
    """
    phases = []
    for time, fit in zip(times, fits):
        phases.append(fit.phase(time))
    return phases

def joined_phase(times, phases, join_time, join_data):
    """ Joins a list of phases, creating a single continuous numpy array.
        This is done by dropping values in the overlap ranges, such that 
        each range has its leftmost overlap dropped.
    """
    pt, pv, tt, tv = peaks.peaks_from_polyfit(join_time, join_data)
    finished=False
    last_tt = 0
    out_time = []
    out_phase = []
    count = 0
    for time, phase in zip(times, phases):
        lidx = 0
        ridx = len(time)
        #trim_t = np.copy(time)
        #trim_ph = np.copy(phase)
        if np.any(time<last_tt):  
            lidx = np.fabs(time-last_tt).argmin()
        try:
            trough_here = np.any(time>=tt[count])
            assert trough_here, "There was a range with no join_data minimum."
        except IndexError:
            finished = True 
        if not finished:
            ridx = np.fabs(time-tt[count]).argmin()
            last_tt = tt[count]
            count = count + 1
        #out_time.append(time[lidx:ridx])
        out_phase.append(phase[lidx:ridx])
    #time_array = np.array(out_time)
    phase_array = np.concatenate(out_phase, axis=0)
    return phase_array
        
def offset_phases(times, fits, join_time, join_data):
    """ Calculate the phases of the power series, then join them together so 
        they agree at the minima of join_data(join_time).  
    """
    pt, pv, tt, tv = peaks.peaks_from_polyfit(join_time, join_data)
    phase_list = phases(times, fits)
    this_min = 0. 
    count = 0
    offset = 0#phase_list[0][0]
    offset_phases = []
    fits_offset = fits[0].phase(tt[0])
    last_tt = 0
    finished=False
    for time, phase, fit in zip(times, phase_list, fits):
        if np.any(time<last_tt):
            offset -= fit.phase(last_tt)
        next_val = phase+offset
        offset_phases.append(next_val)
        try:
            trough_here = np.any(time>=tt[count])
            assert trough_here, "There was a range with no join_data minimum."
        except IndexError:
            finished = True 
        if not finished:
            last_tt = tt[count]
            offset += fit.phase(last_tt)
            count = count+1
    return offset_phases

def divvy_coef(time, data, coef=0.6):
    """Finds the indices d_k at which data is closest to minimal, with the 0th and last
       indices being 0 and len(data). Computes M_k = (d_k + d_{k+1} / 2), 
       W_k = d_{k+1} - d_k, and returns an array such that each row is an index range:
       [[       0          ,  M_0 + coef * W_0],
        [M_1 - coef * W_1  ,  M_1 + coef * W_1],
           .
           .
           .
        [M_N - coef * W_N  ,        -1        ]]  
        The terms 'M_0 -/+ coef * W_0' are rounded to the nearest integer.
        (they must be good indices)
    """
    # Get a sorted list of all extremal indices.
    pk_idx, pk_v, tr_idx, tr_v = peaks.turning_points(data)
    both_idx = [0] + tr_idx + [len(time)]
    both_idx.sort()
    M_k = [(d0 + d1)/2. for d0, d1 in zip(both_idx[:-1], both_idx[1:])]
    fW_k = [coef*(d1 - d0) for d0, d1 in zip(both_idx[:-1], both_idx[1:])]
    
    left = [M - f for M, f in zip(M_k, fW_k)] 
    right = [M + f for M, f in zip(M_k, fW_k)] 
    index_ranges = np.column_stack((left, right))
    index_ranges[0, 0] = 0 
    index_ranges[-1, 1] = len(time)
    return np.around(index_ranges, decimals=0).astype(int)

def divvy(time, data, N_overlap=20):
    """ Breaks 'time' into pieces, each containing about N cycles of
        data, such that each piece overlaps its successor by N_overlap points.
        Returns: a 2D numpy array 'pieces' such that pieces[i, 0] is the start
        and pieces[i, 1] the end index of each piece.
    """
    # Get a sorted list of all extremal indices.
    pk_idx, pk_v, tr_idx, tr_v = peaks.turning_points(data)
    both_idx = pk_idx + tr_idx
    both_idx.sort()
    assert both_idx[-2]+N_overlap < len(time), "failed size assumptions!"
    assert both_idx[-1] < len(time), "failed size assumptions!"
    pieces = np.zeros((len(both_idx)+1, 2))
    pieces[1:, 0] = both_idx[:]
    pieces[:-1, 1] = np.copy(both_idx[:]) + N_overlap
    pieces[-1, 1] = len(time)
    return pieces.astype(int)
