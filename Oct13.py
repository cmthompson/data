# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:29:02 2015

@author: chris
"""

import ramanTools.RamanSpectrum as rs
import UVVistools as uv

def Oct13():
    a= loadtxt('/home/chris/Dropbox/DataWeiss/151013/NMR samples.csv', delimiter = ',', unpack = True, skiprows=1)
    
    for i in range(1,6):
        r = polyfit(a[0,0:100],a[i,0:100],1)
        a[i]-=rs.polyeval(r,a[0])
        peak = uv.findpeak(a[0],a[i],(410,420))
        print peak 
        diameter = -0.000000066521*peak[0]**3+0.00019557*peak[0]**2-0.092352*peak[0]+13.29
        print 'diam', diameter
        epsilon = 21536*diameter**2.3
        print 'eps', epsilon
        print 'conc', 5*peak[1]/epsilon 
        plot(a[0],a[i])
    legend(['10.5','9.5','9.3','8.6','6.8'])
        
    return 0
    
def Oct17():
    a= loadtxt('/home/chris/Dropbox/DataWeiss/151017/thioldots.csv', delimiter = ',', unpack = True, skiprows=1, usecols = (0,6,3,1,2,4,5))
        
    for i in range(1,7):
        r = polyfit(a[0,0:100],a[i,0:100],1)
        a[i]-=rs.polyeval(r,a[0])
        peak = uv.findpeak(a[0],a[i],(410,420))
       #print peak 
        diameter = -0.000000066521*peak[0]**3+0.00019557*peak[0]**2-0.092352*peak[0]+13.29
        #print 'diam', diameter
        epsilon = 21536*diameter**2.3
       # print 'eps', epsilon
        print 'CONC', 5*peak[1]/epsilon 
        plot(a[0],a[i])
    legend(['4.2','6.0',',6.4','7.7','9.7','11.5'])
    
    print 'THIOL EXCHANGE'  
    a= loadtxt('/home/chris/Dropbox/DataWeiss/151017/thiolcapped doits.csv', delimiter = ',', unpack = True, skiprows=1, usecols = (0,1))
    r = polyfit(a[0,0:100],a[1,0:100],1)
    a[1]-=rs.polyeval(r,a[0])
    peak = uv.findpeak(a[0],a[1],(410,420))
   #print peak 
    diameter = -0.000000066521*peak[0]**3+0.00019557*peak[0]**2-0.092352*peak[0]+13.29
    #print 'diam', diameter
    epsilon = 21536*diameter**2.3
   # print 'eps', epsilon
    print 'CONC', 5*peak[1]/epsilon
    
    a= loadtxt('/home/chris/Dropbox/DataWeiss/151017/thioldots.csv', delimiter = ',', unpack = True, skiprows=1, usecols = (0,3,1,2,4,5,6))
        
    for i in range(1,7):
        r = polyfit(a[0,0:100],a[i,0:100],1)
        a[i]-=rs.polyeval(r,a[0])
        peak = uv.findpeak(a[0],a[i],(410,420))
       #print peak 
        diameter = -0.000000066521*peak[0]**3+0.00019557*peak[0]**2-0.092352*peak[0]+13.29
        #print 'diam', diameter
        epsilon = 21536*diameter**2.3
       # print 'eps', epsilon
        print 'CONC', 5*peak[1]/epsilon 
        plot(a[0],a[i])
    legend(['200','400',',800','2000','3200','4000'])
        
    return 0
    
def Oct17NMRfitting():
#    a = loadtxt('/home/chris/Dropbox/DataWeiss/151020/MPAexchange on HCN_1000eq.csv',skiprows = 1, usecols = (0,1), delimiter = ',', unpack = True)
#   `r = RamanSpectrum(pandas.Series(a[1],a[0]))
#    w=1E-5
#    g = [0.05,0.03,0.05,.03,.1,.03,.05,.03,.05,1.38,1.39,1.395,1.405,1.41,1.42,1.425,1.43,1.44,w,w,w,w,w,w,w,w,w,0,0]
#    s = fitspectrum(r,(1.34,1.46), 'xGaussian', g)
#    clf()
#    r.plot()
#    for i in s.peaks: plot(s.x,i)
#    plot(s.x,s.y)
#    print s.areas
#    xlim(1.34,1.49)
#    ylim(-0.01,0.1)
    
 
#    a = loadtxt('/home/chris/Dropbox/DataWeiss/151020/MPAexchange on HCN_100eq.csv',skiprows = 1, usecols = (0,1), delimiter = ',', unpack = True)
#    r = RamanSpectrum(pandas.Series(a[1],a[0]))
#    r.name = ''
#    w=1E-5
#    a = 0.04
#    g = [a,a,a,a,a,a,a,a,a,a,a,1.38,1.385,1.392,1.398,1.405,1.41,1.412,1.42,1.43,1.435,1.44,w,w,w,w,w,w,w,w,w,w,w,0,0]
#   # s =  fitspectrum(r,(1.34,1.46), 'xGaussian', g)
#    
#    r.plot()
#    for i in s.peaks: plot(s.x,i)
#    plot(s.x,s.y)
#    print s.areas
#    xlim(1.34,1.49)
#    ylim(-0.01,0.1)
    
    a = loadtxt('/home/chris/Dropbox/DataWeiss/151020/MPAexchange on HCN_100eq.csv',skiprows = 1, usecols = (0,1), delimiter = ',', unpack = True)
    r = RamanSpectrum(pandas.Series(a[1],a[0]))
    r.name = ''
    w=1E-5
    a = 0.04
    g = [a,0.05,a,a,a,a,a,a,a,a,a,1.725,1.73,1.738,1.743,1.745,1.75,1.755,1.765,1.77,1.78,1.785,w,w,w,w,w,w,w,w,w,w,w,0,0]
    s =  fitspectrum(r,(1.70,1.80), 'xGaussian', g)
    clf()
    r.plot()
    for i in s.peaks: plot(s.x,i)
    plot(s.x,s.y)
    print s.params[0]
    xlim(1.70,1.80)
    ylim(-0.01,0.1)
    return s.areas

def Oct29PPATitrationUVVis():
    a = loadtxt('/home/chris/Dropbox/DataWeiss/151029/151029UVVis.csv', delimiter = ',', unpack=True, skiprows=1)
    for i in range(1,8):
        a[i]-=a[i,0]
        peak = uv.findpeak(a[0],a[i],(410,420))
        print peak
        plot(a[0],a[i])
        diameter = -0.000000066521*peak[0]**3+0.00019557*peak[0]**2-0.092352*peak[0]+13.29
        print 'diam', diameter
        epsilon = 21536*diameter**2.3
       # print 'eps', epsilon
        print 'CONC', 5*peak[1]/epsilon 
    legend(['vial0','vial1','vial2','vial3','vial4','vial5','vial6',])
    ylabel('Absorbance')
    xlabel('wavelength (nm)')
    return 0
    
def Oct30PPATitrationUVVis():
    a = loadtxt('/home/chris/Dropbox/DataWeiss/151030/151030PPAtitration.csv', delimiter = ',', unpack=True, skiprows=1)
    for i in range(1,11):
        a[i]-=a[i,0]
        peak = uv.findpeak(a[0],a[i],(410,420))
        print peak
        plot(a[0],a[i])
        diameter = -0.000000066521*peak[0]**3+0.00019557*peak[0]**2-0.092352*peak[0]+13.29
        #print 'diam', diameter
        epsilon = 21536*diameter**2.3
       # print 'eps', epsilon
        print 'CONC', 5*peak[1]/epsilon 
    legend(['vial0','vial1','vial2','vial3','vial4','vial5','vial6','vial7','vial8','vial9', 'vial10',])
    ylabel('Absorbance')
    xlabel('wavelength (nm)')
    return 0
