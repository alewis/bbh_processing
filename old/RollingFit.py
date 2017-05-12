import numpy as np
import error as err
import cPickle as pickle

#Takes

class RollingFit:
  #Construct from a time series, fed in as separate arguments.
  #The time interval is assumed constant.
  #If a fit window w is supplied, we'll do a fit right away.
  #This class keeps a pickled dict mapping a parameter vectors, window
  #sizes, and 
  def __init__(self, data, times, fitfunc, nparams,fname,w=-1):
    self.t = times
    self.d = data
    self.fitfunc = fitfunc
    self.n = nparams
    self.thefile
    self.fitdict=dict()


  def fit(p,w,n=1):
    from scipy import optimize
    thiskey=(p,w,n)
    if thiskey in self.fitdict:
      print "Key ", thiskey, "exists! Returning it..."
      return self.fitdict[key]
    
    print "Fitting for key ", thiskey, "..."
    
    if self.n != p.shape[0]:
      err("Fit paramater size inconsistency.")
    if w <= 0:
      err("Invalid fit window.")
    for idx in range(0,self.t.shape[0],n):
 
  #PLOTTING

  #Return the plot to some other interface (useful if you want to make a plot
  #of multiple fits)
  def returnplot(p,w,n=1,lab=""):
    import matplotlib.pyplot as plt
    times,fit = self.fitdict(p,w,n)
    if lab=="":
      lab=(p,w,n)
    times,fit = self.fitdict(p,w,n,label=lab)
    return plt.plot(times,fit)
    
  #Show the plot, then save it
  def showplot(path,p,w,n=1,xlab="Time",ylab="Freq",title=""):
    import matplotlib.pyplot as plt
    times,fit = self.fitdict(p,w,n)
    if title=="":
      title=(p,w,n)
    plt.plot(times,fit)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plot.title(title)
    plt.show()
    plt.savefig(path)


