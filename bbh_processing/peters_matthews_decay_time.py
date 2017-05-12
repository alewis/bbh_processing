"""
This integrates the Peters (1964) equations for two-body inspiral to determine
the Omega_phi(t=0) and r_A(t=0) resulting in a merger t=TMERGE_TARGET M with
initial eccentricity ECC_TARGET.
"""
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from stdout_redirected import stdout_redirected

###Parameter targets
Q = 1 #Mass ratio Q = m1/m2
ECC_TARGET = 0.6
T_MERGE_TARGET = 13000.
T_MAX = 30000
TSTEP = 0.01
D0MIN = 10
D0MAX = 70 
EPSILON = 0.0001

###Constants
G = 1.  #Newton's constant
C = 1.  #Speed of light
M = 1.  #Total mass M = m1 + m2
COEF = G**3.*M/C**5. * Q/(Q+1)**2. #Appears in the ODE's
PHI = (1. + 5.**0.5)/2.

def da_dt(ecc, a_sem):
    """
    Secular decay of the semi-major axis.
    """
    result = COEF*(-64./5.)/(a_sem**3.*(1.-ecc**2)**(7./2.))
    result *= (1. + 73./24. * ecc**2. + 37./96. * ecc**4.)
    return result

def de_dt(ecc, a_sem):
    """
    Secular decay of the eccentricity.
    """
    result = -304./15. * ecc * COEF / (a_sem**4.*(1.-ecc**2.)**(5./2.))
    result *= (1.+(121./304.) * ecc**2.)
    return result

def problem_vector(args, time):
    """
    Returns the ODE system [de_dt(ecc,r_p), dr_dt(ecc,r_p)]
    """
    vec = [de_dt(args[0], args[1]), da_dt(args[0], args[1])]
    return vec

def solve_equations(ecc0, a_0, time):
    """
    This solves the ODE system given by problem_vector with initial
    conditions ecc(t=0) = ecc0, a(t=0) = a_0.
    Output: a 2D numpy array (soln[:,0] - e(t), soln[:,1] - a(t).
    """
    y_0 = [ecc0, a_0]
    #with stdout_redirected():
    soln = odeint(problem_vector, y_0, time)
    return soln

def plot_system(ecc0, a_0, time):
    """Solves the equations and plots the results.
    """
    soln = solve_equations(ecc0, a_0, time)
    idx_of_first_zero = np.nonzero(soln[:, 1])[0][-1] + 1
    make_plot(time, soln[:, 1], soln[:, 0], idx_of_first_zero)


def eccentricity(r_a, omega_0):
    """
    Return the Keplerian eccentricity for a given r_A and omega_0
    """
    return 1. - omega_0**2 * r_a**3 / G*M

def apastron_freq(r_a, ecc):
    """
    Return orbital frequency at apastron from eccentricity and r_A.
    """
    return np.sqrt(G*M*(1.-ecc)/r_a**3.)

def compute_tmerge(ecc, r_a0, time, plot=True):
    """
    This finds tmerge as a function of eccentricity and initial separation.
    """
    a_0 = r_a0/(1.+ecc)
    #time = np.arange(0, t_max, tstep)
    soln = solve_equations(ecc, a_0, time)
    ecc_t = soln[:, 0]
    a_t = soln[:, 1]
    idx_of_first_zero = np.nonzero(a_t)[0][-1] + 1


    if len(a_t) - idx_of_first_zero <= 1:
        print "MERGER NOT FOUND"
        return -1 
    else:
        tmerge = time[idx_of_first_zero]
        #print "TMERGE: ", tmerge
    if plot:
        make_plot(time, a_t, ecc_t, idx_of_first_zero)
    return tmerge

def make_plot(t, a_t, ecc, idx_of_first_zero=None):
    """Plot the merger vs. time.
    """
    t = t[:idx_of_first_zero]
    a_t = a_t[:idx_of_first_zero]
    ecc = ecc[:idx_of_first_zero]

    r_a = a_t*(1. + ecc)
    r_p = a_t*(1. - ecc)
    plt.plot(t, r_a, label="r_a", color="green")
    plt.plot(t, a_t, label="a", color="blue")
    plt.plot(t, r_p, label="r_p", color="green")
    ax = plt.gca()
    cur_ylim = ax.get_ylim()
    ax.set_ylim(0, cur_ylim[1])
    ax.set_xlabel("t")
    ax.set_ylabel("r")
    plt.legend(loc="lower left")
    plt.show()

def tmerge_error(ecc, a0, t_target, time):
    """returns T=(t_merge(ecc,d0)-t_target)**2
    """
    t_merge = compute_tmerge(ecc, a0, time, plot=False)
    #sqr_error = np.fabs(t_merge - t_target)
    print "t_merge = ", t_merge
    sqr_error = (t_merge - t_target)**2
    return sqr_error

def golden_section_search(initial_bound, ecc, t_target, t_max, tstep):
    """Does golden section search on tmerge_error to find optimal d0.
    """
    thisa = initial_bound[0]
    thisb = initial_bound[1]
    time = np.arange(0., t_max, tstep)
    done = False
    while done==False:
        c_1 = (PHI - 1.)*thisa + (2. - PHI)*thisb
        c_2 = (2. - PHI)*thisa + (PHI - 1.)*thisb
        f_c1 = np.sqrt(tmerge_error(ecc, c_1, t_target, time))
        f_c2 = np.sqrt(tmerge_error(ecc, c_2, t_target, time))
        print "c1= ", c_1, "c2= ", c_2
        print "f(c1)= ", f_c1, "f(c2)= ", f_c2
        if f_c1 == f_c2:
            print "Break condition: f_c1==f_c2."
            root = c_1
            done = True
        if c_2 - c_1 < EPSILON:
            print "Break condition: inside EPSILON"
            if f_c1 < f_c2:
                root = c_1
                done = True
            else:
                root = c_2
                done = True
        if f_c1 < f_c2:
            thisb = c_2
        else:
            thisa = c_1
    compute_tmerge(ecc, root, time)
    return root 
def find_initial_params(d0min=10, d0max=100, t_merge_target=12500,
                        ecc_target=0.4, tstep=0.1):
    """Find omega_0 and d0 giving ecc = ecc_target and t_m = t_merge_target.
       Uses golden section search to minimize
      T=(t_merge(ecc,d0)-t_merge_target)**2, which has a unique minimum at
      t_merge=t_merge_target (and is continuous).
    """
    t_right = T_MAX
    #a0min = d0min/(1+ecc_target)
    #a0max = d0max/(1+ecc_target)
    initial_bound = (d0min, d0max)
    d0_opt = golden_section_search(initial_bound, ecc_target, t_merge_target,
                                   t_right, tstep)
    #d0_opt = a0_opt*(1+ecc_target)
    print "d0 = ", d0_opt
    omega_opt = apastron_freq(d0_opt, ecc_target)
    print "OMEGA0 = ", omega_opt

if __name__ == "__main__":
    #compute_tmerge(ECC_TARGET, 100, 100000, 0.01)
    find_initial_params(D0MIN, D0MAX, T_MERGE_TARGET, ECC_TARGET, TSTEP)
