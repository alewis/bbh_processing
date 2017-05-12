#!/usr/bin/env python
#coding: utf-8

""" Filtered and non-filtered Wigner Ville Distribution 
    This work is licensed under a Creative Commons Attribution 3.0 Unported License.
    Frank Zalkow, 2013

    ported from tfrwv from the matlab/octave Time-Frequency Toolbox
    (http://tftb.nongnu.org/)
    differences to tfrwv:
       - only takes one signal -> only the auto wigner distribution,
         no cross wigner distribution
       - argument make_analytic for calculating the wigner distribut-
         ion of the analytic signal (less crossterms because negative
         frequency content is zero) """

from time import clock
from numpy import abs, arange, shape, array, ceil, zeros, conj, ix_, transpose, append, fft, real, float64, linspace, sqrt
from scipy.signal import hilbert
from scipy import interpolate
from math import log, ceil, floor

""" get next power of 2 """
def nextpow2(n):
    return 2 ** ceil(log(n, 2))

""" auxilliary function for wvd if trace is true """
def disprog(i, N, steps):
    global begin_time_disprog 
    if i == 0:
        begin_time_disprog = clock()
    if i == (N-1):
        print "100 %% complete in %f seconds." % (clock() - begin_time_disprog)
        del begin_time_disprog
    elif (floor(i * steps / float(N)) != floor((i-1) * steps / float(N))) :
        print "%d" % (floor(i * steps / float(N)) * ceil(100.0 / float(steps))),

""" calculate the wigner ville distribution of a signal"""
def wvd(signal, t=None, N=None, trace=0, make_analytic=True):
    if make_analytic:
        x = hilbert(signal)
    else:
        x = array(signal)
        
    if x.ndim == 1: [xrow, xcol] = shape(array([x]))
    else: raise ValueError("Signal x must be one-dimensional.")
        
    if t is None: t = arange(len(x))
    if N is None: N = len(x)
    
    if (N <= 0 ): raise ValueError("Number of Frequency bins N must be greater than zero.")
    
    if t.ndim == 1: [trow, tcol] = shape(array([t]))
    else: raise ValueError("Time indices t must be one-dimensional.")
    
    if xrow != 1:
        raise ValueError("Signal x must have one row.")
    elif trow != 1:
        raise ValueError("Time indices t must have one row.")
    elif nextpow2(N) != N:
        print "For a faster computation, number of Frequency bins N should be a power of two."
    
    tfr = zeros([N, tcol], dtype='complex')
    if trace: print "Wigner-Ville distribution",
    for icol in xrange(0, tcol):
        ti = t[icol]
        #if ti==0:
        #  raise ValueError("ti==0 causes problems.")
        taumax = min([ti, xcol-ti-1, int(round(N/2.0))-1])
        tau = arange(-taumax, taumax+1)
        indices = [int((N+ttau)%N) for ttau in tau]
        linds = [int(ti+ttau) for ttau in tau]
        rinds = [int(ti-ttau) for ttau in tau]
        
        lhsinds = ix_(indices,[icol])
        tfr[lhsinds] = transpose(array(x[linds] * conj(x[rinds]), ndmin=2))
        tau=int(round(N/2))+1
        if ((ti+1) <= (xcol-tau)) and ((ti+1) >= (tau+1)):
            if(tau >= tfr.shape[0]): tfr = append(tfr, zeros([1, tcol]), axis=0)
            tfr[ix_([tau], [icol])] = array(0.5 * (x[ti+tau] * conj(x[ti-tau]) + x[ti-tau] * conj(x[ti+tau])))
        if trace: disprog(icol, tcol, 10)
    
    tfr = real(fft.fft(tfr, axis=0))
    f = 0.5*arange(N)/float(N)
    return (transpose(tfr), t, f )

""" get the filtered wvd by multiplying the wvd and the stft """ 
def filtered_wvd(wvd, stft):
    qstft = abs(stft)
    qstft = float64(qstft * qstft)    
    
    bigstft = zeros(shape(wvd[0]), float64)
    
    x = arange(0, shape(qstft)[0])
    y = arange(0, shape(qstft)[1])
    
    xx = linspace(x.min(), x.max(), shape(wvd[0])[0])
    yy = linspace(y.min(), y.max(), shape(wvd[0])[1])
    
    interpolator = interpolate.RectBivariateSpline(x,y,qstft, kx=1,ky=1)
    
    bigstft = interpolator(xx,yy)
    
    return (sqrt(abs(bigstft * wvd[0])), wvd[1], wvd[2])
