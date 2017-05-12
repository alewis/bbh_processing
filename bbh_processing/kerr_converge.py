from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import bbh_processing.analyze_kerr_data as kerr
import bbh_processing.equatorial_freq_methods as eqfreq
import cPickle as pickle

prs=OptionParser()
prs.add_option("--infile", type="string", help="Name of kerr orbit.")
prs.add_option("--outfile", type="string", help="Where to pickle things.")
prs.add_option("--eq", action="store_true", default=False)

(opts, args) = prs.parse_args()

def save_ommag(infile, outfile): 
    time, ommag = kerr.even_sampled_omega(infile, ptsperorbit=1000, method="savgol")
    outpath = outfile+"ommag.p"  
    pickle.dump([time, ommag], open(outpath, "wb"))    
    print "Saved ommag"

def save_omr(time, ommag, Ns, outfile):
    omrs = [eqfreq.omega_r_sliding(time, ommag, N, refine=False, minidx=10) for N in Ns]
    outpath = outfile+"omr.p"
    pickle.dump([Ns, omrs], open(outpath, "wb"))
    print "Saved omr"

def save_omph(time, ommag, Ns, outfile):
    omphs = [eqfreq.omega_ph_sliding(time, ommag, N, refine=False, minidx=10) for N in Ns]
    outpath = outfile+"omph.p"
    pickle.dump([Ns, omphs], open(outpath, "wb"))
    print "Saved omph"

def save_omth(time, costheta, ommag, Ns):
    omths = [unwrap.omega_theta(time, costheta, time, ommag, N) for N in Ns]
    outpath = outfile+"omth.p"
    pickle.dump([Ns, omths], open(outpath, "wb"))
    print "Saved omth"

def save_omrho(infile, outfile, Ns):
    rraw = kerr.loadfile(infile)
    time = rraw[:, 0]
    rvec = rraw[:, 1:]
    rvec[:, 3] = 0.
    ommag_rho = savgol_omegafromr(time, rvec, 7, delta=(time[1]-time[0]))
    omphs = [eqfreq.omega_ph_sliding(time, ommag_rho, N, refine=False, minidx=10) for N in Ns]
    outpath = outfile+"omrho.p"
    pickle.dump([Ns, omphs, ommag_rho], open(outpath, "wb"))
    print "Saved omrho"
    
def do_converge(infile, outfile, eq):
    save_ommag(infile, outfile)
    time, ommag = pickle.load(open(outfile+"ommag.p", "rb"))
    
    Ns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    save_omr(time, ommag, Ns, outfile)
    save_omph(time, ommag, Ns, outfile)

    if not eq:
        rvec = kerr.loadfile(infile)[:, 1:]
        costheta = kerr.costheta(rvec, [0, 0, 1])
        save_omth(time, costheta, ommag, Ns)

    save_omrho(infile, outfile, Ns)


if __name__=="__main__":
    do_converge(opts.infile, opts.outfile, opts.eq)
