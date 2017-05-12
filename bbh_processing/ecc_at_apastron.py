#Generate a list of eccentricities at the apastrons of a data set,
#and store them along with the associated orbital frequencies in a file
def make_ecc_dict(
                   sepname,
                   pname,
                   mode    = "proper"
                  ):
 
  #Load in orbital frequency data 
  f = open(pname)
  Omega,OmegaMag=pickle.load(f)
  f.close()
  cols=()
  
  #Load in proper or coordinate separation depending on choice of mode 
  if mode=="proper":
    cols = (0,1)
  elif mode=="coord":
    cols = (0,2)
  else:
    raise ValueError("Mode must be 'coord' or 'proper'")
  t, r = np.loadtxt(sepname,usecols=cols,unpack=True)
  
  #trim the noisy bit off the beginning
  t = t[200:]
  r = r[200:]
  turns, turn_vals, rngs = turning_points(r,first_left=0,first_right=1,stride=1)
  
  #trim so we start at apastron
  if (turn_vals[0]<turn_vals[1]):
    turns = turns[1:]
    turn_vals = turn_vals[1:]
    rngs = rngs[1:]
  
  #and finish at periastron
  if (turn_vals[-2]<turn_vals[-1]):
    turns = turns[:-2]
    turn_vals=turn_vals[-2]
    rngs=rngs[-2]
  ecclist = []
  rlist = []
  freqlist = []
  for i in range(0,len(turns),2):
    ecclist.append(ecc(r[turns[i]],r[turns[i+1]]))
    print OmegaMag[turns[i]]
    freqlist.append(OmegaMag[:,1][turns[i]])
    rlist.append(r[turns[i]]) 
  print rlist, freqlist 
  #zip r and freq into duples to use as keys 
  keys = zip(rlist,freqlist)

  #and return the dictionary
  return dict(zip(keys,ecclist))
