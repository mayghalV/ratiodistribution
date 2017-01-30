#Motivation: A regession has predicted a value for variable X (mu_x) with error e_x
#Similarly, another independent regression has predicted a value for value Y (mu_y) with error e_y

#X = mu_x +- e_x
#Y = mu_y +- e_y

#Question: If we define Z:=X/Y. What is the error bounds of Z? ie e_z?

#See: https://en.wikipedia.org/wiki/Ratio_distribution#Gaussian_ratio_distribution
#for mathematical information

import numpy as np

from scipy.integrate import quad as _quad
from scipy.stats import norm as _norm
from scipy.optimize import fsolve as _fsolve

import matplotlib.pyplot as _plt

from math import sqrt as _sqrt
from math import exp as _exp
from math import pi as _pi

#=========================================================  a
def _a(sd_x, sd_y, z):
    return _sqrt(z**2/(sd_x**2)+1/(sd_y**2))

#=========================================================  b
def _b(mu_x, sd_x, mu_y, sd_y, z):
    return mu_x/(sd_x**2)*z+mu_y/(sd_y**2)

#=========================================================  c
def _c(mu_x, sd_x, mu_y, sd_y):
    return (mu_x**2)/(sd_x**2)+(mu_y**2)/(sd_y**2)

#=========================================================  d
def _d(mu_x, sd_x, mu_y, sd_y, z):
    a_val = _a(sd_x, sd_y, z)
    b_val = _b(mu_x, sd_x, mu_y, sd_y,z)
    c_val = _c(mu_x, sd_x, mu_y, sd_y)
    return _exp((b_val**2-c_val*a_val**2)/(2*(a_val**2)))

#=========================================================  sdFromError
def sdFromError(error, confidence):
    #Standard deviation from error bound and confidence interval (assumes normal distribution)
    ci=(1-confidence)/2
    return error/_norm.ppf(1-ci)

#=========================================================  pdf
def pdf(z, mu_x, sd_x, mu_y, sd_y):
    
    a_val = _a(sd_x, sd_y, z)                  #a,b,c,d
    b_val = _b(mu_x, sd_x, mu_y, sd_y,z)
    c_val = _c(mu_x, sd_x, mu_y, sd_y)
    d_val = _d(mu_x, sd_x, mu_y, sd_y,z)
    
    part1 = b_val*d_val/((a_val**3)*_sqrt(2*_pi)*sd_x*sd_y)
    part1 *= _norm.cdf(b_val/a_val)-_norm.cdf(-b_val/a_val)
    
    part2 = 1/(a_val**2*_pi*sd_x*sd_y)*_exp(-c_val/2)
    
    return part1+part2

#=========================================================  cdf
def cdf(z, mu_x, sd_x, mu_y, sd_y):
    #Find P(Z<=z)
    return _quad(pdf, -np.inf, z, args=(mu_x, sd_x, mu_y, sd_y))[0]

#=========================================================  inverse cdf
def invCdf(cdf_val, mu_x, sd_x, mu_y, sd_y):
    #Find z such that P(Z<=z)=cdf_val
    
    #Inner function shifting to zero
    def newP_z(Z, mu_x, sd_x, mu_y, sd_y):
        return cdf(Z, mu_x, sd_x, mu_y, sd_y)-cdf_val
    
    return _fsolve(newP_z, mu_x/mu_y, args = (mu_x, sd_x, mu_y, sd_y))[0]

#=============================================================== error bounds of Z
def error(mu_x, e_x, mu_y, e_y, confidence=0.95):
    #Confidence interval 2 sided (if 95%, ci = (100%-95%)/2)
    #Returns 95% confidence error bounds
    
    #Calculations
    sd_x = sdFromError(e_x,confidence)                   #Calculate standard deviations from error bounds
    sd_y = sdFromError(e_y,confidence)
    ci=(1-confidence)/2
    
    
    exp_z = mu_x/mu_y       #Expected value of z
    e_u = invCdf(1-ci,mu_x, sd_x, mu_y, sd_y)-exp_z       #Error upper
    e_l = exp_z - invCdf(ci, mu_x, sd_x, mu_y, sd_y)      #Error lower
    
    mid = str(exp_z)
    
    print('Bounds: [' + mid + '-' + str(e_l) + ',' + mid + '+' + str(e_u) +']')
    
    return [exp_z-e_l, e_u+exp_z]
#=============================================================== plot pdf
def plotpdf(mu_x, e_x, mu_y, e_y, confidence=0.95):
    
    sd_x = sdFromError(e_x,confidence)   #Calculate standard deviations from error bounds
    sd_y = sdFromError(e_y,confidence)
    
    #Plot pdf
    x_val = np.linspace(mu_x/mu_y-max(sd_x,sd_y),mu_x/mu_y+max(sd_x,sd_y),100)
    y_val = [pdf(x,mu_x, sd_x, mu_y, sd_y) for x in x_val]
    
    _plt.plot(x_val, y_val)
    _plt.show()
    return
    