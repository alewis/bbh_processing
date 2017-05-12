#Return N values of q evenly spaced in symmetric mass ratio
#from q=qmin to q=qmax.
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
def eta(q):
    return q/(q+1.)**2.
def symm_mass(eta_min, eta_max, N):
    q_0s = np.linspace(20., 1., 10000)
    eta_q = eta(q_0s) 
    q_eta = UnivariateSpline(eta_q, q_0s, k=3, s=0)
    plt.plot(eta_q, q_0s)
    plt.show()

    # eta_min = eta(qmin)  
    # eta_max = eta(qmax)
    etas = np.linspace(eta_min, eta_max, N)
    qs = q_eta(etas)
    print etas/0.25, qs 

if __name__ == "__main__":
    symm_mass(0.25, 0.25/4., 5.)
    
    

