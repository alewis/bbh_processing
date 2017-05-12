#plot_omega_component.py
#Creates a .dat file plotting time against the frequency
#of oscillations in the r,phi,or theta direction.

import argparse

parser = argparse.ArgumentParser(description="Using a 'Horizons.h5'"+
   "and 'Omega.pkl' file storing apparent horizon and orbital"+ 
    "frequency data from a bbh run, create .dat files storing the"+
    "component frequencies (r, theta, phi) as functions of time."+
    "Each goes in a separate file.")

parser.add_argument('-r',action='store_true',default=False,
    help="Create 'r' output if flagged.")
parser.add_argument('-theta',action='store_true',default=False,
    help="Create 'theta' output if flagged.")
parser.add_argument('-phi',action='store_true',default=False,
    help="Create '-phi' output if flagged.")
parser.add_argument('--hpath',default="Horizons.h5",
    help="Path to the Horizons.h5 file.") 
parser.add_argument('--ompath', default="Omega.pkl",
    help="Path to the Omega.pkl file.") 
parser.add_argument('--outpath', default="omega",
    help="Tag added to the output's filename (full name is"+
          "{outpath}_{component}.dat.") 

from bbh_processing.combined_tools import \
    omega_r, omega_theta, omega_phi
import numpy as np

def make_omega_theta(hname,outname):
  #plot omega_theta 
  t, om = omega_theta(hname)
  
  out=file(outname,'w')
  out.write("""#
    # [1] = t
    # [2] = omega_theta 
    """)
  np.savetxt(out,np.transpose([t,om]))

def make_omega_phi(om_name,outname):
  #plot omega_phi
  t, om = omega_phi(om_name)
  out=file(outname,'w')
  out.write("""#
    # [1] = t
    # [2] = omega_phi 
    """)
  np.savetxt(out,np.transpose([t,om]))

def make_omega_r(om_name,outname):
  #plot omega_r 
  t, om = omega_r(om_name)
  
  out=file(outname,'w')
  out.write("""#
    # [1] = t
    # [2] = Omega_r 
    """)
  np.savetxt(out,np.transpose([t,om]))


args = parser.parse_args()
if (not args.r and not args.theta and not args.phi):
  print "No output flag specified; program will do nothing."

if (args.r): 
  rpath = args.outpath+"_r.dat"
  print "Making omega_r file '",rpath,"'..."
  make_omega_r(args.ompath,rpath)
  print "\tMade. \n"
else:
  print "-r not specified; skipping r..."

if (args.theta): 
  thpath = args.outpath+"_th.dat"
  print "Making omega_theta file '",thpath,"'..."
  make_omega_theta(args.hpath,thpath)
  print "\tMade. \n"
else:
  print "-theta not specified; skipping theta..."

if (args.phi): 
  phipath = args.outpath+"_phi.dat"
  print "Making omega_phi file '",phipath,"'..."
  make_omega_phi(args.ompath,phipath)
  print "\tMade. \n"
else:
  print "-phi not specified; skipping phi..."
