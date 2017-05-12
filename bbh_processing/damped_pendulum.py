#d^2(th)/dt^2 + u d(th)/dt + w^2 sin(th) = 0
import numpy as np
from scipy import integrate 
import matplotlib.pyplot as plt



def spring(mu, x0, gamma, y, t):
    xdot = y[1]
    x = y[0]
    xdotdot = -mu*xdot*np.cos(x)*((x0-np.sin(x))**gamma) - np.sin(x)
    
    #*((1.0/ - x
    #xhat = x0 + x 
    #term1 = xdot#*((1.0/(xhat)**1.) 
    #term2 = x 
    #xdotdot = -term1 - term2 
    return [xdot, xdotdot]

def damped_pendulum(v_ini, x_ini, tf, mu=0.2, x0=1.5, gamma=-2.):
    a_t = np.arange(0, tf, 0.001)
    thisspring = lambda Y, t: spring(mu, x0, gamma, Y, t)
    asol = integrate.odeint(thisspring, [x_ini, v_ini], a_t)
    astack = np.c_[a_t, asol[:, 0], asol[:, 1]]
    return astack[:, 0], astack[:, 1], astack[:, 2]


def gw_pend_work(mu, x_s, y, t):
    th = y[0]
    thdot = y[1]
    thdotdot = -mu * thdot * np.cos(th)*((x_s - np.sin(th))**5.)
    thdotdot = thdotdot - np.sin(th)
    return [thdot, thdotdot]

def gw_pendulum(v_ini, x_ini, tf, mu=0.2, x_s=1.5): 
    a_t = np.linspace(0, tf, 10000)
    thispend = lambda Y, t: gw_pend_work(mu, x_s, Y, t)
    asol = integrate.odeint(thispend, [x_ini, v_ini], a_t)
    return a_t, asol[:, 0], asol[:, 1]
