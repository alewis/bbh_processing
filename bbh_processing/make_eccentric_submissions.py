import itertools
from subprocess import call
from sys import exit
from string import replace
import os
import re

MAKE_ID = False
ACTUALLY_RUN = False

OLDMINLEV = "$MinLev=1;"
NEWMINLEV = "$MinLev=3;"

OLDMAXLEV = "$MaxLev=3;"
NEWMAXLEV = "$MaxLev=3;"

#Whether to send emails to the SpEC database (generally you don't want this).
#0 for no, 1 for yes.
OLDTRACKRUN = "$TrackRunInDatabase = 1;"
NEWTRACKRUN = "$TrackRunInDatabase = 0;"

#The tags for the input dictionary
ECC="ecc"
Q="q"
OMEGA0="Omega0"
ADOT="adot"
D0="D0"
INC="inc"
CHIA="chiA"
CHIAX="chiAx"
CHIAY="chiAy"
CHIAZ="chiAz"
CHIB="chiB"
CHIBX="chiBx"
CHIBY="chiBy"
CHIBZ="chiBz"

#Just a wrapper around the Python 'call' function, which runs bash commands.
def cmd(the_command):
  call(the_command,shell=True)

#Reads the input file into a dictionary. Involves a little but of ugly string
#manipulation. 
def read_input_file(name="eccentric_submissions.input"):
  with open(name,"r") as the_file:
    lines = the_file.readlines()
  #str_list = [l.split()[2:][0].split(",") for l in lines]
  #split along = with whitespace and newlines removed 
  str_list = [l.replace(" ", "").strip().split('=') for l in lines]

  in_dict = {} 
  for key,l in str_list:
    floatl = [float(s) for s in l.split(",")] 
    in_dict[key]=floatl 
  in_dict_error_trapping(in_dict)
  return in_dict


def in_dict_error_trapping(in_dict):
  qsz = len(in_dict[Q])

  eccsz = len(in_dict[ECC])
  orblists = [in_dict[OMEGA0],in_dict[ADOT],in_dict[D0]]
  if any(eccsz!=len(x) for x in orblists): 
    raise ValueError("Orbital parameters lists had different sizes.")

  incsz = len(in_dict[INC])
  spinlists = [in_dict[CHIAX], in_dict[CHIAY],in_dict[CHIAZ]]
  spinlists += [in_dict[CHIBX], in_dict[CHIBY],in_dict[CHIBZ]]
  if any(incsz!=len(x) for x in spinlists): 
    raise ValueError("Spin parameter lists had different sizes.")

  sz=str(qsz*eccsz*incsz) 
  answer = raw_input("This will submit " + sz + " runs; 'y' to confirm: ") 
  if answer != 'y':
    print "Input was not 'y' - aborting..."
    exit()


#Return a string containing the contents of a Params.input file
def params_file_string(
    Omega0 = 0.0,
    adot0 = 0.0,
    D0 = 0.0,
    MassRatio = 1,
    SpinA = (0,0,0),
    SpinB = (0,0,0),
    Evolve = 1,
    UseSKS = 1):
  the_string = "# Params file created by make_eccentric_submissions.py\n" 
  the_string += "# Set the initial data parameters\n"
  the_string += "\n"
  
  the_string += "# Orbital parameters\n"
  the_string += "$Omega0 = " + str(Omega0) + ";\n"
  the_string += "$adot0 = " + str(adot0) + ";\n"
  the_string += "$D0 = " + str(D0) + ";\n"
  the_string += "\n"
  
  the_string += "# Physical parameters\n"
  the_string += "$MassRatio = " + str(MassRatio) + ";\n"
  the_string += "@SpinA = " + "(" + str(SpinA[0]) + "," + str(SpinA[1]) + ","
  the_string += str(SpinA[2]) + ");\n"
  the_string += "@SpinB = " + "(" + str(SpinB[0]) + "," + str(SpinB[1]) + ","
  the_string += str(SpinB[2]) + ");\n"
  the_string += "\n"
  
  the_string += "# Evolve after initial data completes?\n"
  the_string += "$Evolve = " + str(Evolve) + ";\n"
  the_string += "# Use SKS flag? (if false, uses CF)\n"
  the_string += "$UseSKS = " + str(UseSKS) + ";\n"
  return the_string

def make_orbd(orb):
  orbd={}
  orbd[OMEGA0]=orb[1]
  orbd[ADOT]=orb[2]
  orbd[D0]=orb[3]
  return orbd
def make_incd(inc):
  print inc
  incd={}
  incd[CHIA]=inc[1]
  incd[CHIB]=inc[2]
  return incd

#Generates a list of tuples (dirname,param_string), where dirname is the 
#top-level name of the directory where the simulation will run and 
#param_string is the contents of the appropriate Params.input file.
def generate_dirnames_and_params(in_dict):
  qs = in_dict[Q]
  orbzp = zip(in_dict[ECC],in_dict[OMEGA0],in_dict[ADOT],in_dict[D0])
  
  spnzp = zip(in_dict[INC],
              zip(in_dict[CHIAX],in_dict[CHIAY],in_dict[CHIAZ]),
              zip(in_dict[CHIBX],in_dict[CHIBY],in_dict[CHIBZ]))
  f_dict={}
  strlist = [] 
  rep = lambda x : str(x).replace(".","_")
  for q,orb,inc in itertools.product(qs,orbzp,spnzp):
    qstr = rep(q)
    eccstr = rep(orb[0])
    incstr = rep(inc[0])
    s = "q"+qstr+"_ecc"+eccstr+"_inc"+incstr
    orbd = make_orbd(orb)
    incd = make_incd(inc)
    param_string = params_file_string(
        Omega0 = orbd[OMEGA0],
        adot0 = orbd[ADOT],
        D0 = orbd[D0],
        MassRatio = q,
        SpinA = incd[CHIA],
        SpinB = incd[CHIB],
        Evolve = 1,
        UseSKS = 1)
    strlist.append((s,param_string)) 
  return strlist 

def directory_loop(function,arg,dirnames,basepath="./"):
  paths = [basepath + d for d in dirnames]
  for path in paths:
    os.chdir(path)
    function(arg)
  os.chdir(basepath)

#Creates a list of dirs with names 'dirnames' from path 'basepath', 
#then runs PrepareID inside each.
#If specpath is supplied, will use -d specpath as the SpEC repo. 
#If reduce_ecc is true, will run with eccentricity reduction.
def make_dirs_and_prepare_ID(dirnames,basepath="./",specpath=None,reduce_ecc=False,
    verbose=True):
  specopt = ""
  if specpath is not None:
    specopt = "-d " + specpath  
  
  if reduce_ecc:
    eccopt = " -reduce-ecc "
  else:
    eccopt = " -no-reduce-ecc "

  print "Creating directories..."
  mkdirs = ["mkdir " + basepath + d for d in dirnames]
  for mkdir in mkdirs: 
    if verbose: print mkdir
    cmd(mkdir) 

  print "Preparing ID..."
  command = "PrepareID -t bbh" + eccopt + specopt
  directory_loop(cmd,"touch B",dirnames,basepath)
  #cds = ["cd " + basepath + d for d in dirnames]
  #for cd in cds:
  #  cmd(cds)
    #cmd(command)

#Input: three lists, one of qs, one of eccs, one of incs. Makes strings
#"qA_eccB_incC", one for each permutation.
def fname_strings_q_ecc_inc(qlist,ecclist,inclist):
  strlist = [] 
  rep = lambda x : str(x).replace(".","_")
  for q,ecc,inc in itertools.product(qlist,ecclist,inclist):
    qstr = rep(q) 
    eccstr = rep(ecc)
    incstr = rep(inc)
    s = "q"+qstr+"_ecc"+eccstr+"_inc"+incstr
    strlist.append(s)
  return strlist

#Actually create the Params.input files.
def write_param_files(dlist,plist,basepath="./"):
  print "Writing Param files..."
  for d,par in zip(dlist,plist):
    with open(basepath+d+"/Params.input","w") as the_file:
      the_file.write(par)

#Edit Ev/DoMultipleRuns.input
def edit_do_multiple_runs(dlist,basepath):
  paths = basepath+dlist+"/Ev/DoMultipleRuns.input"
  for p in paths: 
    with open(p,"r") as source:
      data = source.read() 
    new_data=replace(data,OLDMINLEV,NEWMINLEV)
    new_data=replace(data,OLDMAXLEV,NEWMAXLEV)
    new_data=replace(data,OLDTRACKRUN,NEWTRACKRUN)
    with open(p,"w") as dest:
      dest.write(new_data)

#Submit the runs.
def submit(dlist,basepath):
  directory_loop(cmd,"qsub ./Submit.sh",dlist,basepath)

if __name__=="__main__":
  thisdir = os.getcwd() 
  basepath = thisdir+"/test/"
  in_dict=read_input_file()
  dlist,plist = generate_dirnames_and_params(in_dict) 
  
  make_dirs_and_prepare_ID(dlist,basepath)
  write_param_files(dlist,plist,basepath)
  edit_dot_multiple_runs(dlist,basepath)
  if ACTUALLY_RUN:
    submit(dlist,basepath) 
  else:
    print "ACTUALLY_RUN was false; will not submit."
  

