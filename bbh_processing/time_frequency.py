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
import numpy as np
from numpy.lib import stride_tricks
import matplotlib.pylab as plt
from scipy.signal import hilbert
from scipy import interpolate
from math import log, ceil, floor
def plotstft(data, binsize=1000, plotpath=None, colormap="jet",
    samplerate=2):
    s = stft(data, binsize)
    
    sshow, freq = logscale_spec(s, factor=1.0, sr=samplerate)
    ims = 20.*np.log10(np.abs(sshow)/10e-6) # amplitude to decibel
    
    timebins, freqbins = np.shape(ims)
    
    plt.figure(figsize=(15, 7.5))
    plt.imshow(np.transpose(ims), origin="lower", aspect="auto", cmap=colormap, interpolation="none")
    plt.colorbar()

    plt.xlabel("time (s)")
    plt.ylabel("frequency (hz)")
    plt.xlim([0, timebins-1])
    plt.ylim([0, freqbins])

    xlocs = np.float32(np.linspace(0, timebins-1, 5))
    plt.xticks(xlocs, ["%.02f" % l for l in ((xlocs*len(data)/timebins)+(0.5*binsize))/samplerate])
    ylocs = np.int16(np.round(np.linspace(0, freqbins-1, 10)))
    plt.yticks(ylocs, ["%.02f" % freq[i] for i in ylocs])
    
    if plotpath:
        plt.savefig(plotpath, bbox_inches="tight")
    else:
        plt.show()
        
    plt.clf()

def logscale_spec(spec, sr=2, factor=20.):
    timebins, freqbins = np.shape(spec)

    scale = np.linspace(0, 1, freqbins) ** factor
    scale *= (freqbins-1)/max(scale)
    scale = np.unique(np.round(scale))
    
    # create spectrogram with new freq bins
    newspec = np.complex128(np.zeros([timebins, len(scale)]))
    for i in range(0, len(scale)):
        if i == len(scale)-1:
            newspec[:,i] = np.sum(spec[:,scale[i]:], axis=1)
        else:        
            newspec[:,i] = np.sum(spec[:,scale[i]:scale[i+1]], axis=1)
    
    # list center freq of bins
    allfreqs = np.abs(np.fft.fftfreq(freqbins*2, 1./sr)[:freqbins+1])
    freqs = []
    for i in range(0, len(scale)):
        if i == len(scale)-1:
            freqs += [np.mean(allfreqs[scale[i]:])]
        else:
            freqs += [np.mean(allfreqs[scale[i]:scale[i+1]])]
    
    return newspec, freqs

""" short time fourier transform of audio signal """
def stft(sig, frameSize, overlapFac=0.5, window=np.hanning):
    win = window(frameSize)
    hopSize = int(frameSize - np.floor(overlapFac * frameSize))
    
    # zeros at beginning (thus center of 1st window should be for sample nr. 0)
    samples = np.append(np.zeros(np.floor(frameSize/2.0)), sig)    
    # cols for windowing
    cols = np.ceil( (len(samples) - frameSize) / float(hopSize)) + 1
    # zeros at end (thus samples can be fully covered by frames)
    samples = np.append(samples, np.zeros(frameSize))
    
    frames = stride_tricks.as_strided(
        samples, shape=(cols, frameSize), 
        strides=(samples.strides[0]*hopSize, samples.strides[0])).copy()
    frames *= win
    
    return np.fft.rfft(frames)   
"""get next power of 2 """
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

""" calculate the wigner ville distribution of a matrix """
def wvd(data, t=None, N=None, trace=0, make_analytic=True):
    if make_analytic:
        x = hilbert(data[1])
    else:
        x = array(data[1])
        
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
        raise ValueError("Time indicies t must have one row.")
    elif nextpow2(N) != N:
        print "For a faster computation, number of Frequency bins N should be a power of two."
    
    tfr = zeros([N, tcol], dtype='complex')
    if trace: print "Wigner-Ville distribution",
    for icol in xrange(0, tcol):
        ti = t[icol]
        taumax = min([ti, xcol-ti-1, int(round(N/2.0))-1])
        tau = arange(-taumax, taumax+1)
        indices = ((N+tau)%N)
        tfr[ix_(indices, [icol])] = transpose(array(x[ti+tau] * conj(x[ti-tau]), ndmin=2))
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
