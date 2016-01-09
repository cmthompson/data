# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:53:24 2015

@author: chris
"""
from scipy.optimize import curve_fit
from ramanTools.RamanSpectrum import *
from numpy import *
from matplotlib.pyplot import *
from copy import deepcopy
import copy

from scipy.optimize import minimize
from matplotlib import gridspec
from UVVistools import findpeak


def Aug27abs():
    figure()
    title('Aug27')
    a = loadtxt('/home/chris/Dropbox/DataWeiss/150827/CdS_synthAug27.csv', skiprows = 1, unpack= True,delimiter=',')
    a[1]-=a[1][0]
    
    plot(a[0],a[1])
    
    b = loadtxt('/home/chris/Dropbox/DataWeiss/150827/CdSsynthAug27_fluor.csv',skiprows = 1, unpack= True,delimiter=',')
    b = RamanSpectrum(pandas.Series(b[3],b[0]))
    b.set_name('myname')
    b[:]/=b[408]
    b.plot()
    b[:]/=b[408]
    r = fitspectrum(b,(407,470),'xGaussian',[1,0.5,0.2,423,439,450,18,10,18,0,0])
    print r.params
    for i in r.peaks:
        plot(r.x,i)
    b.plot()
    print r.params[0][6]
    print 'FWHM fluorescence =', 2*numpy.sqrt(r.params[0][6]*log(2))
    
    figure()
    title('UVVis absorption of CdS quantum dots')
    c = loadtxt('/home/chris/Dropbox/DataWeiss/150828/CdS_synthAug28.csv', skiprows = 1, unpack= True,delimiter=',')
    c[1]-=c[1][0]
    plot(c[0],c[1])
    ylabel('absorbance (a. u.)')
    xlabel('wavelength (nm)')
    annotate('$\lambda_{max}$ = 412 nm \n diameter = 3.8 nm',(412,0.51))
    
    figure()
    title('Aug31')
    d = loadtxt('/home/chris/Dropbox/DataWeiss/150831/CdS_synthAug31.csv', skiprows = 1, unpack= True,delimiter=',')
    d[1]-=d[1][0]
    d[1]/=0.497
    plot(d[0],d[1])
    
    dfluor = loadtxt('/home/chris/Dropbox/DataWeiss/150831/CdS_synthAug31_fluor.csv',skiprows = 1, unpack= True,delimiter=',')
    dfluor = RamanSpectrum(pandas.Series(dfluor[3],dfluor[0]))
    dfluor.set_name('myname')
    dfluor[:]-=min(dfluor)
    dfluor[:]/=dfluor[422]
    dfluor.plot()
    legend(['absoption', 'fluorescence'])
    ylabel('absorbance (a. u.)')
    xlabel('wavelength (nm)')
    
    r = fitspectrum(dfluor,(407,670),'xGaussian',[1,1,423,580,18,200,0,0])
    print r.params
#    for i in r.peaks:
#        plot(r.x,i)

    
    print 'FWHM fluorescence  =', 2*numpy.sqrt(r.params[0][4]*log(2))
 
    
    
    figure()
    title('all spectra corrected for dilution')
    subplot(121)
    a[1]*=(1652/52)
    c[1]*=(1675/73)
    d[1]*=(1731/63)
    
    plot(a[0],a[1])
    plot(c[0],c[1])
    plot(d[0],d[1])
    legend(['Aug27', 'Aug28', 'Aug31'])
    ylabel('absorbance of original solution (a. u.)')
    
    average = (40*a[1]+40*c[1]+80*d[1])/160
    
    subplot(122)
    
    plot(a[0],average)
    legend(['average'])
    return 0
    

          
    
            