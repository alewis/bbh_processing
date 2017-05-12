import multiprocessing as mproc
import numpy as np
import itertools
def make_vector(i,L):
  print i, L
  a = np.empty(L)
  a.fill(i)
  print a.transpose()

def modify_shared_array(i):
  y[:,i] = i

def start_proc():
  print "Starting: ", mproc.current_process().name

def main():
  if __name__ == '__main__':
    procids = list(range(10))
    arg_two = 5
    print "Inputs: ", procids, "Arg two: " , arg_two


    #map(starred_call, itertools.izip(procids,itertools.repeat(arg_two)))

    pool_size = mproc.cpu_count() 
    print pool_size

    pool = mproc.Pool(processes = pool_size,
                      initializer = start_proc,
                      maxtasksperchild=2,
                      )
    pool.map(modify_shared_array,procids)
    pool.close()
    pool.join()
    print y

if __name__=="__main__":
  main()
