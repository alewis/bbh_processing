""" Stores classes representing power series fits to bbh data.
"""
import numpy as np
import matplotlib.pyplot as plt
import warnings
import scipy.optimize as opt
from bbh_processing.unwrap import unwrap

# class KeplerPowerSeries(object):
    # """Represent a function of form
       # z = K * sin(phi)/(1 + e cos(phi) ) + B
       # Where phi is a power series of order o.
    # """
    # def __init__(self, time, data, order=3,
                 # ecc=0.1):
        # self.order = order
        # self.offset = (np.amin(data) + np.amax(data) / 2.)

        # self.c = time[0] + 0.5*(time[-1] - time[0])
        # K, ecc, coefs = self.__fit(time, data, order, ecc, K)
        # self.ecc = ecc
        # self.coefs = coefs
        # self.K = K

    # def __call__(self, time, ecc=None, coefs=None):
        # """Return the full function:
           # f(t) = self.amp * cos(self.phase) + self.offset
        # """
        # if ecc is None:
            # ecc = self.ecc
        # phi = self.power_series(time, coefs)
        # output = self.K * np.sin(phi)/(1. + ecc * np.cos(phi)) + self.offset
        # output = self.K * np.cos(dotted) + self.offset
        # return output

    # def phase(self, time, coefs=None):
        # """Return the 'phi' function.
        # """
        # if coefs is None:
            # p_coefs = self.coefs
        # else:
            # p_coefs = np.poly1d(coefs)
        # phase = p_coefs(time - self.c)
        # return phase

    # def __fit(self, time, data, order, guess=None):
        # err = lambda p, x, y: self(x, p) - y
        # warnings.simplefilter('ignore', np.RankWarning)
        # shift_t = time - self.c
        # if guess is None:
            # period_guess = np.pi/(time[-1] - time[0])
            # phase_guess = np.arccos(self.__arccos_argument(data))
            # guess_coefs = np.polyfit(shift_t, phase_guess, order)
            # poly = np.poly1d(guess_coefs)
        # else:
            # poly = guess.coefs
        # p0, success = leastsq(err, poly.c, args=(time, data))
        # return np.poly1d(p0)

class CosinePowerSeries(object):
    """Represents an amplitude times the cosine of a power series
    plus an offset:
    f(t) = A * cos(a_0 + a_1(t-T) + a_2(t-T)^2) + B
    (is the phase shift T important?)
    """
    
    def __init__(self, time, data, order=3, 
                 method="leastsq", amp=None):
        assert time.shape==data.shape, "Time, data must have same shape."
        self.order = order
        self.c = time[0] + 0.5*(time[-1] - time[0])
        self.__set_amp_offset(data)
        coefs = np.poly1d(self.__fit(time, data, order, method=method))
        self.coefs = coefs

    def __set_amp_offset(self, data):
        """Calculate the amplitude and offset of the function.
        """
        self.amp = (np.amax(data) - np.amin(data)) / 2.
        min_off = np.amin(data)+self.amp
        max_off = np.amax(data)-self.amp
        self.offset = (min_off + max_off) / 2.
        
        return

    def __call__(self, time, coefs=None):
        """Return the full function:
           f(t) = self.amp * cos(self.phase) + self.offset
        """
        phase = self.phase(time, coefs=coefs)
        output = self.amp * np.cos(phase) + self.offset

        return output
    
    def phase(self, time, coefs=None):
        """Return the 'phase' interior to the cosine:
          phi(t) = a_0 + a_1(t) + a_2(t)^2 ...
        """
        if coefs is None:
            p_coefs = self.coefs
        else:
            p_coefs = np.poly1d(coefs)
        phase = p_coefs(time - self.c)
        return phase
    
    def vary_everything(self, time, guesses):
        """Return the full function:
           f(t) = self.amp * cos(self.phase) + self.offset
        """
        guess_copy = list(guesses)
        offset = guess_copy.pop()
        amp = guess_copy.pop()
        phase = self.phase(time, coefs=guess_copy)
        output = amp * np.cos(phase) + offset
        return output
    
    def jac(self, time, coefs):
        """Jacobian of the full function:
          df/da_i = -self.amp*sin(self.phase)*(t-C)^i
        """
        my_phase = self.phase(time, coefs=coefs)
        #the_jac = np.zeros((len(time), len(coefs)+1))
        the_jac = []
        for i in np.arange(len(coefs)):
            the_jac.append(-self.amp*np.sin(my_phase)*((time-self.c)**i)) 
            #the_jac[:, i] = -self.amp*np.sin(my_phase[:])*(time[:]-self.c)**i        
        return the_jac[::-1]

    

    # def unwrapped_phase(self, time, coefs=None):
        # """
          # N is the number of half-cycles into the fit this series represents.
          # We must add Pi * N to the final result to get an unwrapped phase.
        # """
        # phase = self.phase(time, coefs)
        # is_odd = lambda x: x & 0x1
        # if not is_odd(self.N):
            # phase = 2*phase[0] - phase
            # #phase = phase
        # N_orbits = (self.N-(self.N_offset/2))/2
        # unwrapped = phase + N_orbits * 2*np.pi
        # return unwrapped

    # def phase(self, time, coefs=None):
        # """return the 'phase' interior to the cosine:
          # phi(t) = a_0 + a_1(t) + a_2(t)^2 ...
        # """
        # if coefs is None:
            # coefs = self.coefs
        # shift = time # - self.c
        # answer = np.zeros(len(time))
        # for i in range(0, len(time)):
            # this_answer = 0.
            # for j in range(0, len(coefs)):
               # this_answer += coefs[j] * time[i]**self.powers[j]
            # answer[i] = this_answer
        # return answer



    def __str__(self):
        outstring = "Order: " + str(self.order) + "\n"
        outstring += "Coefs: " + str(self.coefs) + "\n"
        outstring += "C: " + str(self.c) + "\n"
        outstring += "Amp: " + str(self.amp) + "\n"
        outstring += "Offset: " + str(self.offset) + "\n"
        return outstring

    def __arccos_argument(self, data):
        """ Invert the functional form and restrict to (-1, 1) so that arccos
        gives a finite answer.
        """
        arg = (data - self.offset) / self.amp
        too_big = arg > 1.
        arg[too_big] = 1.
        too_small = arg < -1.
        arg[too_small] = -1.
        return arg

    def __cost_function(self, coefs, time, data):
        """ Return sum(t) [data - f(t)]^2.
        """
        residuals = data - self(time, coefs=coefs)
        cost = np.sum(np.square(residuals))
        #print "Cost:",cost
        return cost

    def __cost_gradient(self, coefs, time, data):
        """ Return the gradient of the cost function w.r.t the fit params.
            This is grad_i = -2*sum_t(resid * df/da_i).
            df/da_i = -A * sin[phase(t-C)] * (t - C) ^i.
            so grad_i = 2*A * sum_t(resid * (t-C)^i * sin[phase(t-C)])
        """
        #print "Coefs: ", coefs
        #print "Order: ", self.order
        #print time
        residuals = data - self(time, coefs=coefs)
        sin_terms = np.sin(self.phase(time, coefs=coefs))
        product = residuals*sin_terms

        shift_t = time-self.c
        sums = np.zeros(self.order+1)
        for i in np.arange(self.order+1):
            sums[i] = np.sum(product * shift_t**i)
        derivs = 2*self.amp*sums
            
       # print "Grad: ", derivs
        #print "Total: ", np.sqrt(np.sum(np.square(derivs)))
        return derivs

    def __cost_hessian(self, coefs, time, data):
        """ Return the Hessian of the cost function w.r.t the fit params.
            This is H_ij = -2 * sum_t(df/da_j + df/da_i da_j).
            df/da_j = -A * sin[phase(t-C)] * (t - C) ^ j.
            df/da_i da_j = -A * cos[phase(t-C)] * (t-C)^i * (t-C)^j.
            so H_ij = 2*A*sum_t(
                        (t-C)^j * (sin[phase(t-C)] + cos[phase(t-C)]*(t-C)^i)
                        )
            'Data' is not actually necessary, but is required by scipy to
            be one of the arguments.
        """
        phase = self.phase(time, coefs=coefs)
        sin_terms = np.sin(phase)
        cos_terms = np.cos(phase)
        shift_t = time-self.c
        hessian = np.zeros((self.order+1, self.order+1))
        for i in range(self.order+1):
            for j in range(self.order+1):
                pwr_j = shift_t**j
                pwr_i = shift_t**i
                these_summands = pwr_j * (sin_terms + pwr_i * cos_terms)
                hessian[i, j] = 2*self.amp*np.sum(these_summands)
        return hessian
    
    def __make_guesses(self, time, data, order):
        """Generate coef. guesses by taking the arccos of the data.
        """
        shift_t = time - self.c
        phase = np.arccos(self.__arccos_argument(data))
        unwrap_outs = unwrap(phase, junk_inds=0)
        phase_guess = unwrap_outs[0]
        flipped = unwrap_outs[1]
        if not ((len(phase_guess)==len(flipped) and len(shift_t)==len(flipped))):
            raise ValueError("Internal size mismatch! Array shapes: "+
                  "\n\tShift_t: "+str(shift_t.shape)+
                  "\n\tFlipped: "+str(flipped.shape)+
                  "\n\tPhase_guess: "+str(phase_guess.shape))
        guess_coefs = np.polyfit(shift_t, phase_guess, order)
        poly = np.poly1d(guess_coefs)
        # plt.plot(time, phase)
        # plt.plot(time, phase_guess)
        # plt.show()
        return poly

    def __newton_fit(self, time, data, guess_poly):
        """Fit using the Newton-CG method. This requires analytic gradient
           and hessian.
           THIS METHOD SEEMS NOT TO WORK AT PRESENT, SINCE THE GRADIENT FAILS
        """
        cost_check = lambda x: self.__cost_function(x, time, data)
        grad_check = lambda x: self.__cost_gradient(x, time, data)
        #check = opt.check_grad(cost_check,
        #                                      grad_check,
        #                                      guess_poly.c
        #                                     )
        result = opt.minimize(cost_check,
                               guess_poly.c,
                               method="Newton-CG",
                               jac=grad_check
                               #hess=self.__cost_hessian
                              )
        # result = opt.minimize(self.__cost_function,
                              # guess_poly.c,
                              # args=(time, data),
                              # method="Newton-CG",
                              # jac=self.__cost_gradient
                              # #hess=self.__cost_hessian
                             # )
        coefs = result.x
        # check = opt.check_grad(cost_check,
                                              # grad_check,
                                              # coefs
                                             # )
        if not result.success:
            print "WARNING: optimize failed!"
            print "Error message: ", result.message
        return coefs

    def __fit(self, time, data, order, guess=None, method="leastsq"):
        """Do a fit to the power series coefficients. Currently, leastsq and
           newton-CG fit methods are implemented; the latter takes advantage
           of the fact that the gradient and hessian of f w.r.t to the fit
           params. can be computed analytically.
        """
        warnings.simplefilter('ignore', np.RankWarning)
        if guess is None:
            poly = self.__make_guesses(time, data, order)
        else:
            poly = guess.coefs
        
        if method == "newton":
            coefs = self.__newton_fit(time, data, poly)
        elif method == "leastsq":
            err = lambda p: self(time, p) - data
            jaclam = lambda p: self.jac(time, p)
            #jaclam = None
            result = opt.leastsq(err, poly.c, Dfun=jaclam, col_deriv=True)
            coefs = result[0]
        elif method == "fit_everything":
            my_guess = list(poly.c) + [self.amp, self.offset]
            err = lambda p: self.vary_everything(time, p) - data
            result = list(opt.leastsq(err, my_guess)[0])
            self.offset = result.pop()
            self.amp = result.pop()
            coefs = result
        elif method == "arccos":
            coefs = poly.c
        else:
            raise ValueError("Invalid method!")
        return np.poly1d(coefs)
