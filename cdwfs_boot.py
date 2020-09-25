#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 09:27:42 2020

@author: yan
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 08:58:04 2020

@author: yan
"""

import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits

from scipy.integrate import quad
from scipy.optimize import curve_fit

import matplotlib as mpl
from matplotlib.pyplot import plot

#########################################

ener = np.arange(0.075,15.0,0.15)


y_lir_nox=[]
area_lir_nox=[]
width_lir_nox=[]

p1_lir_nox=[]
reda=[]

def func2(x, a, b):
    return a*x+b




f=fits.open('f_irh.fits')

index=np.array([6,7,8,9,10,11,12,13,14,15,16])

xdata = f[0].data[0][32:55]
new_x = np.delete(xdata, index)


f=fits.open('f_irh_boot.fits')


def func(x,c, d,e):
    para = popt_pl[0]
    parb = popt_pl[1]

    return c*np.exp(-0.5*((x-e)/d)**2)+para*x+parb

def integrand(x):
    return popt[0]*np.exp(-0.5*((x-popt[2])/popt[1])**2)

p1=np.array([0.5,0.3])
p2=np.array([0.5,0.3,6.4])


width=[]

for i in range(0,500):
    
    y =f[0].data[i][32:55]
    ydata = y #+ y_noise


    new_y = np.delete(ydata, index)

    popt_pl, pcov = curve_fit(func2, new_x, new_y, p0=p1)

    popt, pcov = curve_fit(func, xdata, ydata, p0=p2,
                           bounds=((0, 0,6.4), (np.inf, 0.5,6.7)))


    fungau = popt[0]*np.exp(-0.5*((xdata-popt[2])/(popt[1]))**2)
    funpl = popt_pl[0]*xdata+popt_pl[1]
    funmod = fungau+funpl


    ans, err = quad(integrand, 5., 8.)
    y_line=popt_pl[0]*popt[2]+popt_pl[1]
    y_lir_nox.append(y_line)
    area_lir_nox.append(ans)
    width_lir_nox.append(ans/y_line)    
    #print popt
    width.append(width_lir_nox)

    plot(xdata, fungau, "b--", linewidth=1)
    plot(xdata, funmod, "-", color="purple",linewidth=4,label="Best Fits")
    plot(xdata, y, 'b.',label="Net Flux")
    
    i=i+1



#plot.set_ylim(0.01,1.25)
#plot.set_xlim(5.0, 7.99)

