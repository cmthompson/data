# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 13:05:06 2015

@author: chris
"""

import weissdatavariables
from scipy.optimize import curve_fit
from ramanTools.RamanSpectrum import *
from numpy import *
from matplotlib.pyplot import *
from UVVistools import indivQY,findpeak, HWHM
from copy import deepcopy
import copy
from scipy.optimize import minimize
from matplotlib import gridspec
import os



    
def anthracenespectra():
     figure()
     ax1 = subplot(121)
     ax2 = subplot(122)
     
     a = loadtxt('150928/150928UVVis.csv',delimiter = ',', unpack = True, skiprows = 1,usecols=(0,8))
     b = loadtxt('150929/150929UVVis.csv',delimiter = ',', unpack = True, skiprows = 1,usecols=(0,8))
     c = loadtxt('150930/150930UVVis.csv',delimiter = ',', unpack = True, skiprows = 1,usecols=(0,13))
     d = loadtxt('151001/151001UVVis.csv', delimiter = ',', unpack = True, skiprows = 1,usecols=(0,7))
     e = loadtxt('151002/151002UVVis.csv', delimiter = ',', unpack = True, skiprows = 1,usecols=(0,20))
     f = loadtxt('151003/151003UVVis.csv', delimiter = ',', unpack = True, skiprows = 1,usecols=(0,19))
     g = loadtxt('151005/151005UVVis.csv', delimiter = ',', unpack = True, skiprows = 1,usecols=(0,19))
     a = RamanSpectrum(pandas.Series(a[1][::-1],a[0][::-1]))
     b = RamanSpectrum(pandas.Series(b[1][::-1],b[0][::-1]))
     c = RamanSpectrum(pandas.Series(c[1][::-1],c[0][::-1]))
     d = RamanSpectrum(pandas.Series(d[1][::-1],d[0][::-1]))
     e = RamanSpectrum(pandas.Series(e[1][::-1],e[0][::-1]))
     f = RamanSpectrum(pandas.Series(f[1][::-1],f[0][::-1]))
     g = RamanSpectrum(pandas.Series(g[1][::-1],g[0][::-1]))
     
     
     afluor = RamanSpectrum(pandas.Series(*loadtxt('150928/150928fluor/anthracene.dat',unpack = True, skiprows = 1,usecols=(3,0))))
     bfluor = RamanSpectrum(pandas.Series(*loadtxt('150929/150929fluor/anthracene.dat',unpack = True, skiprows = 1,usecols=(3,0))))
     cfluor = RamanSpectrum(pandas.Series(*loadtxt('150930/150930fluor/anthracene.dat',unpack = True, skiprows = 1,usecols=(3,0))))
     dfluor = RamanSpectrum(pandas.Series(*loadtxt('151001/151001fluor/trial1/anthracene.dat',unpack = True, skiprows = 1,usecols=(3,0))))
     efluor = RamanSpectrum(pandas.Series(*loadtxt('151002/151002fluor/anthracene.dat',unpack = True, skiprows = 1,usecols=(3,0))))
     ffluor = RamanSpectrum(pandas.Series(*loadtxt('151003/151003fluor/anthracene.dat',unpack = True, skiprows = 1,usecols=(3,0))))
     gfluor = RamanSpectrum(pandas.Series(*loadtxt('151005/151005fluor/anthracene.dat',unpack = True, skiprows = 1,usecols=(3,0))))
     
     absorbance = array([])
     fluorescence= array([])
     for i in (a,b,c,d,e,f,g,):
        # i.autobaseline((314,390),order= 0)
         #i.plot(ax=ax1)
         i.smoothbaseline((290,300),(390,400),_plot=False,ax=ax1)
         
         print i[350]/i[400]
         print i[350]
         i.plot(ax=ax1)
         absorbance=numpy.append(absorbance,i[350])
     ax1.legend(['0','1','2','3','4','5','7'])
     ax1.set_ylim(0,0.1)
     
     
     fluorlist = (afluor,bfluor,cfluor,dfluor,efluor,ffluor,gfluor,)
#     for i in range(len(fluorlist)):
#         fluorlist[i][:]/=absorbance[i]
     fluorescence= array(list((i[420]*78.2032212661 for i in fluorlist)))
     for i in fluorlist:
        # fluorescence=numpy.append(fluorescence,i[420]*78.2032212661)
         i.plot(ax=ax2)
     ax2.legend(['0','1','2','3','4'])
     figure()
     plot((1-10**(-absorbance))/fluorescence)
     
     
     return 0

   
    
def trial3():
    """"Decrease in Quantum yield of Quantum dots CdS capped iwth PPA over time.  Trial 3"""
 
    sample1=list()
    sample2=list()
    sample3=list()
    sample4=list()
    sample5=list()
    sample6=list()
    PPAbuff=list()
    
    figure()
    Aax1 = subplot(231)
    Aax2 = subplot(232)
    Aax3 = subplot(233)
    Aax4 = subplot(234)
    Aax5 = subplot(235)
    Aax6 = subplot(236)
    
    figure()
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(234)
    ax5 = subplot(235)
    ax6 = subplot(236)
   
    #day0 
    sample1.append(indivQY('151001/151001UVVis.csv', 'b','h','151001/151001fluor/trial3/t3day0N2pH7.dat', '151001/151001fluor/trial3/anthracene.dat',UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7'))
    sample2.append(indivQY('151001/151001UVVis.csv', 'c','h','151001/151001fluor/trial3/t3day0N2pH9.dat', '151001/151001fluor/trial3/anthracene.dat',UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9'))
    sample3.append(indivQY('151001/151001UVVis.csv', 'd','h','151001/151001fluor/trial3/t3day0N2pH11.dat', '151001/151001fluor/trial3/anthracene.dat',UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11'))
    sample4.append(indivQY('151001/151001UVVis.csv', 'e','h','151001/151001fluor/trial3/t3day0airpH7.dat', '151001/151001fluor/trial3/anthracene.dat',UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7'))
    sample5.append(indivQY('151001/151001UVVis.csv', 'f','h','151001/151001fluor/trial3/t3day0airpH9.dat', '151001/151001fluor/trial3/anthracene.dat',UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9'))
    sample6.append(indivQY('151001/151001UVVis.csv', 'g','h','151001/151001fluor/trial3/t3day0airpH11.dat', '151001/151001fluor/trial3/anthracene.dat',UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11'))
    #PPAbuff.append(indivQY('151001/151001UVVis.csv', 'u','h','151001/151001fluor/trial3/t3day0PPAbuff.dat', '151001/151001fluor/trial3/anthracene.dat',UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'PPAbuffer'))
    
    #day1
    anthracenefile = '151002/151002fluor/trial3/anthracene.dat'
    uvvisfile = '151002/151002UVVis.csv'
    sample1.append(indivQY(uvvisfile, 'n','u','151002/151002fluor/trial3/t3d1N2pH7.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'o','u','151002/151002fluor/trial3/t3d1N2pH9.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'p','u','151002/151002fluor/trial3/t3d1N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'q','u','151002/151002fluor/trial3/t3d1airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'r','u','151002/151002fluor/trial3/t3d1airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 's','u','151002/151002fluor/trial3/t3d1airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    #PPAbuff.append(indivQY('151002/PPAbuffdots151002.csv', 'b','c','151001/151001fluor/trial3/t3day0PPAbuff.dat', anthracenefile,UVVisplot = Aax1,fluorplot=ax1,day = 0,label = 'PPAbuffer'))
    ##day2
    day=2
    anthracenefile = '151003/151003fluor/anthracene.dat'
 
    uvvisfile = '151003/151003UVVis.csv'
    fluorfolder = '151003/151003fluor/trial3/'
    sample1.append(indivQY(uvvisfile, 'n','t',fluorfolder+'t3N2pH7.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'o','t',fluorfolder+'t3N2pH9.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'p','t',fluorfolder+'t3N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(0)#indivQY(uvvisfile, 'q','t',fluorfolder+'t3airpH7.dat', anthracenefile,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'q','t',fluorfolder+'t3airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'r','t',fluorfolder+'t3airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    #PPAbuff.append(indivQY(uvvisfile, 's','t',fluorfolder +'PPABuffDots.dat', anthracenefile,UVVisplot = Aax1,fluorplot=ax1,day = 0,label = 'PPAbuffer'))
    
    day=4
    anthracenefile = '151005/151005fluor/anthracene.dat'
    uvvisfile = '151005/151005UVVis.csv'
    fluorfolder = '151005/151005fluor/trial3/'
    sample1.append(indivQY(uvvisfile, 'n','t',fluorfolder+'t3N2pH7.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'o','t',fluorfolder+'t3N2pH9.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'p','t',fluorfolder+'t3N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(0)#indivQY(uvvisfile, 'q','t',fluorfolder+'t3airpH7.dat', anthracenefile,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'q','t',fluorfolder+'t3airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'r','t',fluorfolder+'t3airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
  #  PPAbuff.append(indivQY(uvvisfile, 's','u',fluorfolder +'PPABuffDots.dat', anthracenefile,UVVisplot = Aax1,fluorplot=ax1,day = 0,label = 'PPAbuffer'))
  
  
    for p in (Aax1,Aax2,Aax3,Aax4,Aax5,Aax6,ax1,ax2,ax3,ax4,ax5,ax6): 
        p.legend(['0','1','2','3','4','5'])
    figure()
    days = [0,1,2,4]
    plot(days,sample1,'s-',label = 'N2pH7')
    plot(days,sample2,'s-',label = 'N2pH9')
    plot(days,sample3,'s-',label = 'N2pH11')
    plot(days,sample4,'s-',label = 'airpH7')
    plot(days,sample5,'s-',label = 'airpH9')
    plot(days,sample6,'s-',label = 'airpH11')
    plot(PPAbuff,'s-', label = 'PPAbuffered')
    legend()
    ylabel('band edge QY')
    xlabel('day')
    return (sample1, sample2, sample3, sample4, sample5 ,sample6)
    
def trial1():
    """"Decrease in Quantum yield of Quantum dots CdS capped iwth PPA over time.  Trial 1"""
    sample1=list()
    sample2=list()
    sample3=list()
    sample4=list()
    sample5=list()
    sample6=list()
    
    figure()
    Aax1 = subplot(231)
    Aax2 = subplot(232)
    Aax3 = subplot(233)
    Aax4 = subplot(234)
    Aax5 = subplot(235)
    Aax6 = subplot(236)
    
    figure()
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(234)
    ax5 = subplot(235)
    ax6 = subplot(236)
    ##day0
    anthracenefile = '150928/150928fluor/anthracene.dat'
    sample1.append(indivQY('150928/150928UVVis.csv', 'b','i','150928/150928fluor/day0N2pH7.dat',anthracenefile, fluorescencerange=(405,466),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY('150928/150928UVVis.csv', 'c','i','150928/150928fluor/day0N2pH9.dat',anthracenefile, fluorescencerange=(407,460),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY('150928/150928UVVis.csv', 'd','i','150928/150928fluor/day0N2pH11.dat',anthracenefile,fluorescencerange=(410,470), UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY('150928/150928UVVis.csv', 'e','i','150928/150928fluor/day0airpH7.dat',anthracenefile,fluorescencerange=(405,466), UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY('150928/150928UVVis.csv', 'f','i','150928/150928fluor/day0airpH9.dat',anthracenefile,fluorescencerange=(407,460), UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY('150928/150928UVVis.csv', 'g','i','150928/150928fluor/day0airpH11.dat',anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    ##day1
    anthracenefile = '150929/150929fluor/anthracene.dat'
    uvvisfile = '150929/150929UVVis.csv'
    sample1.append(indivQY(uvvisfile, 'b','i','150929/150929fluor/day1N2pH7.dat',anthracenefile,fluorescencerange=(405,466),   UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'c','i','150929/150929fluor/day1N2pH9.dat',anthracenefile,fluorescencerange=(407,460),  UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','i','150929/150929fluor/day1N2pH11.dat',anthracenefile,fluorescencerange=(410,470), UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','i','150929/150929fluor/day1airpH7.dat',anthracenefile,fluorescencerange=(405,466), UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','i','150929/150929fluor/day1airpH9.dat',anthracenefile,fluorescencerange=(407,460), UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','i','150929/150929fluor/day1airpH11.dat',anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    ##day2
    anthracenefile = '150930/150930fluor/anthracene.dat'
    uvvisfile = '150930/150930UVVis.csv'
    sample1.append(indivQY(uvvisfile, 'b','n','150930/150930fluor/day2N2pH7.dat',anthracenefile,fluorescencerange=(405,466),   UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'c','n','150930/150930fluor/day2N2pH9.dat',anthracenefile, fluorescencerange=(407,460), UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','n','150930/150930fluor/day2N2pH11.dat',anthracenefile,fluorescencerange=(410,470), UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','n','150930/150930fluor/day2airpH7.dat',anthracenefile,fluorescencerange=(405,466), UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','n','150930/150930fluor/day2airpH9.dat',anthracenefile,fluorescencerange=(407,460), UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','n','150930/150930fluor/day2airpH11.dat',anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    
    ##day3
    anthracenefile = '151001/151001fluor/trial1/anthracene.dat'
    uvvisfile = '151001/151001UVVis.csv'
    sample1.append(indivQY(uvvisfile, 'i','h','151001/151001fluor/trial1/t1day3N2pH7.dat',anthracenefile,fluorescencerange=(405,466),  UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'j','h','151001/151001fluor/trial1/t1day3N2pH9.dat',anthracenefile,fluorescencerange=(407,460), UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'k','h','151001/151001fluor/trial1/t1day3N2pH11.dat',anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'l','h','151001/151001fluor/trial1/t1day3airpH7.dat',anthracenefile,fluorescencerange=(405,466), UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'm','h','151001/151001fluor/trial1/t1day3airpH9.dat',anthracenefile,fluorescencerange=(407,460), UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'n','h','151001/151001fluor/trial1/t1day3airpH11.dat',anthracenefile,fluorescencerange=(410,470), UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    ##day4
    anthracenefile = '151002/151002fluor/trial1/anthracene.dat'
 
    uvvisfile = '151002/151002UVVis.csv'
    sample1.append(indivQY(uvvisfile, 'b','u','151002/151002fluor/trial1/t1d4N2pH7.dat',anthracenefile,fluorescencerange=(405,466), fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'c','u','151002/151002fluor/trial1/t1d4N2pH9.dat',anthracenefile,fluorescencerange=(407,460),fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','u','151002/151002fluor/trial1/t1d4N2pH11.dat', anthracenefile,fluorescencerange=(410,470),fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','u','151002/151002fluor/trial1/t1d4airpH7.dat',anthracenefile,fluorescencerange=(405,466),fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','u','151002/151002fluor/trial1/t1d4airpH9.dat', anthracenefile,fluorescencerange=(407,460),fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','u','151002/151002fluor/trial1/t1d4airpH11.dat', anthracenefile,fluorescencerange=(410,470),fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    ##day5
    day=5
    anthracenefile = '151003/151003fluor/anthracene.dat'
 
    uvvisfile = '151003/151003UVVis.csv'
    fluorfolder = '151003/151003fluor/trial1/'
    sample1.append(indivQY(uvvisfile, 'b','t',fluorfolder+'t1N2pH7.dat', anthracenefile,fluorescencerange=(405,466),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'c','t',fluorfolder+'t1N2pH9.dat', anthracenefile,fluorescencerange=(407,460),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','t',fluorfolder+'t1N2pH11.dat', anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','t',fluorfolder+'t1airpH7.dat', anthracenefile,fluorescencerange=(405,466),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','t',fluorfolder+'t1airpH9.dat', anthracenefile,fluorescencerange=(407,460),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','t',fluorfolder+'t1airpH11.dat', anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
     
      ##day7
    day=7
    anthracenefile = '151005/151005fluor/anthracene.dat'
 
    uvvisfile = '151005/151005UVVis.csv'
    fluorfolder = '151005/151005fluor/trial1/'
    sample1.append(indivQY(uvvisfile, 'b','t',fluorfolder+'t1N2pH7.dat', anthracenefile,fluorescencerange=(405,466),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'c','t',fluorfolder+'t1N2pH9.dat', anthracenefile,fluorescencerange=(407,460),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','t',fluorfolder+'t1N2pH11.dat', anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','t',fluorfolder+'t1airpH7.dat', anthracenefile,fluorescencerange=(405,466),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','t',fluorfolder+'t1airpH9.dat', anthracenefile,fluorescencerange=(407,460),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','t',fluorfolder+'t1airpH11.dat', anthracenefile,fluorescencerange=(410,470),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
        
    Aax1.legend(list((str(i) for i in range(6))))
    figure()
    days = [0,1,2,3,4,5,7]
    plot(days,sample1,'bs-',label = 'N2pH7')
    plot(days,sample2,'rs-',label = 'N2pH9')
    plot(days,sample3,'ks-',label = 'N2pH11')
    plot(days,sample4,'bo-',label = 'airpH7')
    plot(days,sample5,'ro-',label = 'airpH9')
    plot(days,sample6,'ko-',label = 'airpH11')
    legend()
    ylabel('band edge QY')
    xlabel('day')
    return (sample1,sample2, sample3,sample4,sample5,sample6)
    
    
def trial2():
    """"Decrease in Quantum yield of Quantum dots CdS capped iwth PPA over time.  Trial 2"""
    ## day 0
    sample1=list()
    sample2=list()
    sample3=list()
    sample4=list()
    sample5=list()
    sample6=list()
    
    figure()
    Aax1 = subplot(231)
    Aax2 = subplot(232)
    Aax3 = subplot(233)
    Aax4 = subplot(234)
    Aax5 = subplot(235)
    Aax6 = subplot(236)
    
    figure()
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(234)
    ax5 = subplot(235)
    ax6 = subplot(236)
    ##day0
    anthracenefile = '150930/150930fluor/anthracene.dat'
    sample1.append(indivQY('150930/150930UVVis.csv', 'h','n','150930/150930fluor/tri2N2pH7.dat',  anthracenefile,UVVisplot = Aax1, fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY('150930/150930UVVis.csv', 'i','n','150930/150930fluor/tri2N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY('150930/150930UVVis.csv', 'j','n','150930/150930fluor/tri2N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY('150930/150930UVVis.csv', 'k','n','150930/150930fluor/tri2airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY('150930/150930UVVis.csv', 'l','n','150930/150930fluor/tri2airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY('150930/150930UVVis.csv', 'm','n','150930/150930fluor/tri2airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    ##day1
    anthracenefile = '151001/151001fluor/trial2/anthracene.dat'
    uvvisfile = '151001/151001UVVis.csv'
    sample1.append(indivQY(uvvisfile, 'o','h','151001/151001fluor/trial2/t2day1N2pH7.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'p','h','151001/151001fluor/trial2/t2day1N2pH9.dat', anthracenefile,UVVisplot = Aax2, fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'q','h','151001/151001fluor/trial2/t2day1N2pH11.dat',anthracenefile,UVVisplot = Aax3, fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'r','h','151001/151001fluor/trial2/t2day1airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 's','h','151001/151001fluor/trial2/t2day1airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 't','h','151001/151001fluor/trial2/t2day1airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
     #day2
    anthracenefile = '151002/151002fluor/trial2/anthracene.dat'
    uvvisfile = '151002/151002UVVis.csv'
    sample1.append(indivQY(uvvisfile, 'h','u','151002/151002fluor/trial2/t2d2N2pH7.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'i','u','151002/151002fluor/trial2/t2d2N2pH9.dat', anthracenefile,UVVisplot = Aax2, fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'j','u','151002/151002fluor/trial2/t2d2N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'k','u','151002/151002fluor/trial2/t2d2airpH7.dat',anthracenefile,UVVisplot = Aax4, fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'l','u','151002/151002fluor/trial2/t2d2airpH9.dat',anthracenefile,UVVisplot = Aax5, fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'm','u','151002/151002fluor/trial2/t2d2airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    ##day3
    day=3
    anthracenefile = '151003/151003fluor/anthracene.dat'
 
    uvvisfile = '151003/151003UVVis.csv'
    fluorfolder = '151003/151003fluor/trial2/'
    sample1.append(indivQY(uvvisfile, 'h','t',fluorfolder+'t2N2pH7.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'i','t',fluorfolder+'t2N2pH9.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'j','t',fluorfolder+'t2N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'k','t',fluorfolder+'t2airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'l','t',fluorfolder+'t2airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'm','t',fluorfolder+'t2airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    ##day3
    day=5
    anthracenefile = '151005/151005fluor/anthracene.dat'
 
    uvvisfile = '151005/151005UVVis.csv'
    fluorfolder = '151005/151005fluor/trial2/'
    sample1.append(indivQY(uvvisfile, 'h','t',fluorfolder+'t2N2pH7.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample2.append(indivQY(uvvisfile, 'i','t',fluorfolder+'t2N2pH9.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'j','t',fluorfolder+'t2N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'k','t',fluorfolder+'t2airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'l','t',fluorfolder+'t2airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'm','t',fluorfolder+'t2airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    Aax1.legend(list((str(i) for i in range(6))))
    figure()
    plot(sample1,'o-',label = 'N2pH7')
    plot(sample2,'o-',label = 'N2pH9')
    plot(sample3,'o-',label = 'N2pH11')
    plot(sample4,'s-',label = 'airpH7')
    plot(sample5,'s-',label = 'airpH9')
    plot(sample6,'s-',label = 'airpH11')
    legend()
    ylabel('band edge QY')
    xlabel('day')
    
    return (sample1, sample2, sample3, sample4, sample5, sample6)
    
    
def brandnewdotsOct6():
    """ fluorescence quantum yields of new CdS dots and the PPA capped version of them"""
    sample1 = list()
    figure()
    ax1 = subplot(211)
    ax2 = subplot(212)
    
    anthracenefile = '151006/151006fluor/anthracene.dat'
    uvvis = '151006/fluorescenceyields.csv'
    print uv.indivQY(uvvis, 'f','i','151006/151006fluor/oleate capped in hexanes.dat','151006/151006fluor/anthracene.dat',
            fluorescencerange = (400,460), excitationwavelength=350,UVVisplot = ax1,fluorplot = ax2,nliq=1.375)
    print uv.indivQY(uvvis, 'd','i','151006/151006fluor/PPA-capped in water.dat','151006/151006fluor/anthracene.dat',
            fluorescencerange = (404,477), excitationwavelength=350,UVVisplot = ax1,fluorplot = ax2,_plot_standard=True)
   
    ax2.legend(['oleate capped', 'PPA-capped'])
    ax2.set_ylabel('differential QY')
    ax1.set_xlabel('wavelength (nm)')
    return 0
    
def brandnewdotsOct7():
    """Fluorescence of oleate, PPA, and MPA capped CdS"""
    sample1 = list()
    figure()
    ax1 = subplot(211)
    ax2 = subplot(212)
    
    sample1.append(indivQY('151006/fluorescenceyields.csv', 'f','i','151006/151006fluor/oleate capped in hexanes.dat',  '151006/151006fluor/anthracene.dat',UVVisplot = ax1,fluorplot = ax2,day = 0,label = 'oleate capped',fluorescencerange = (400,470),nliq=1.375))
    
    anthracenefile = '151007/151007fluor/anthracene.dat'
    sample1.append(indivQY('151007/151007.csv', 'b','e','151007/151007fluor/PPAcappedH2O.dat', anthracenefile,UVVisplot = ax1,fluorplot = ax2,day = 0,fluorescencerange = (406,470)))
    sample1.append(indivQY('151007/151007.csv', 'c','e','151007/151007fluor/PPAcapKOHHCl.dat', anthracenefile,UVVisplot = ax1,fluorplot = ax2,day = 0,label = 'PPAcapped',color = 'r',fluorescencerange = (406,470)))
    sample1.append(indivQY('151007/151007.csv', 'd','e','151007/151007fluor/MPAcappeddots.dat', anthracenefile,UVVisplot = ax1,fluorplot = ax2,day = 0,label = 'MPAcapped',color = 'r',fluorescencerange = (409,455)))
    
    
    for z in range(4):
        i = ax1.lines[z]
        r = findpeak(i.get_xdata(),i.get_ydata(),(410,420))
        print r
        i.set_ydata(i.get_ydata()/r[1])
        if z!=2:
            pass
            #savetxt('home/chris/Desktop/UVVisExchangedDots'+i.get_label()+str(z)+'.csv', transpose([i.get_xdata(),i.get_ydata()]),delimiter = ',')

    ax1.legend(['oleate', 'PPA H2O', 'PPA salt','MPA'])
    ax2.legend(['oleate', 'PPA H2O', 'PPA salt','MPA'])
    ax2.set_ylabel('differential QY')
    ax1.set_xlabel('wavelength (nm)')
    
    print 'PPA fluorescence to MPA fluorescence', sample1[1]/sample1[3]
    figure()
    plot(sample1, 's')
  
    return 0

def Oct12throughx():
    day=5
    
     ## day 0
    sample1=list()
    sample2=list()
    sample3=list()
    sample4=list()
    sample5=list()
    sample6=list()
    
    midgapemission = list()
    mg2=list()
    mg3=list()
    mg4=list()
    
    figure()
    Aax1 = subplot(231)
    Aax2 = subplot(232)
    Aax3 = subplot(233)
    Aax4 = subplot(234)
    Aax5 = subplot(235)
    Aax6 = subplot(236)
    
    figure()
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(234)
    ax5 = subplot(235)
    ax6 = subplot(236)
    
    ax1.set_title('dot1N2')
    ax2.set_title('dot2N2')
    ax3.set_title('dot1air')
    ax4.set_title('dot2air')
    
    halfwidths = list()
    anthracenefile = '151012/151012fluor/anthracene.dat'
    
    uvvisfile = '151012/151012UVVis.csv'
    fluorfolder = '151012/151012fluor'
    sample1.append(indivQY(uvvisfile, 'g','l',fluorfolder+'/dot1N2.dat', anthracenefile,UVVisplot = Aax1, fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = (406,474)))
    sample2.append(indivQY(uvvisfile, 'h','l',fluorfolder+'/dot2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = (406,474)))
    sample3.append(indivQY(uvvisfile, 'i','l',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = (406,474)))
    sample4.append(indivQY(uvvisfile, 'j','l',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = (406,474)))
    sample5.append(indivQY(uvvisfile, 'm','l',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (412,470)))
    
    halfwidths.append(HWHM(loadtxt(uvvisfile,unpack = True, delimiter = ',', skiprows = 1, usecols = (0,6)), (408,420)))
    
    fr = (474,680)
    midgapemission.append(indivQY(uvvisfile, 'g','l',fluorfolder+'/dot1N2.dat', anthracenefile,UVVisplot = Aax1, fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange =fr))
    mg2.append(indivQY(uvvisfile, 'h','l',fluorfolder+'/dot2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = fr))
    mg3.append(indivQY(uvvisfile, 'i','l',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = fr))
    mg4.append(indivQY(uvvisfile, 'j','l',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange =fr))
    #sample5.append(indivQY(uvvisfile, 'm','l',fluorfolder+'dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (412,470)))
  
    
    
    
    ############################################################
    anthracenefile = '151013/151013fluor/anthracene.dat'
    
    uvvisfile = '151013/151013UVVis.csv'
    fluorfolder = '151013/151013fluor/'
    sample1.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dots1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = (406,474)))    
    sample2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dots2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = (406,474)))
    sample3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = (406,474)))
    sample4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = (406,474)))
    sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dot2swithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    halfwidths.append(HWHM(loadtxt(uvvisfile,unpack = True, delimiter = ',', skiprows = 1, usecols = (0,1)), (408,420)))
    fr = (474,680)
    midgapemission.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dots1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = fr))
    mg2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dots2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = fr))
    mg3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = fr))
    mg4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = fr))
    #sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dot2swithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
        
    ###################################################
    anthracenefile = '151014/151014fluor/anthracene.dat'
    
    uvvisfile = '151014/151014UVVis.csv'
    fluorfolder = '151014/151014fluor/'
    sample1.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dot1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = (406,474)))
    sample2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dot2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = (406,474)))
    sample3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = (406,474)))
    sample4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = (406,474)))
    sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    halfwidths.append(HWHM(loadtxt(uvvisfile,unpack = True, delimiter = ',', skiprows = 1, usecols = (0,1)), (408,420)))
    fr = (474,680)
    midgapemission.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dot1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = fr))
    mg2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dot2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = fr))
    mg3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = fr))
    mg4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = fr))
    #sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    
     ###################################################
    anthracenefile = '151015/151015fluor/anthracene.dat'
    
    uvvisfile = '151015/151015UVvis.csv'
    fluorfolder = '151015/151015fluor/'
    sample1.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dot1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = (406,474)))
    sample2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dot2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = (406,474)))
    sample3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = (406,474)))
    sample4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = (406,474)))
    sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    halfwidths.append(HWHM(loadtxt(uvvisfile,unpack = True, delimiter = ',', skiprows = 1, usecols = (0,1)), (408,420)))
    fr = (474,680)
    midgapemission.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dot1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = fr))
    mg2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dot2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = fr))
    mg3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = fr))
    mg4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = fr))
    #sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    
       ###################################################
    anthracenefile = '151016/151016fluor/anthracene.dat'
    
    uvvisfile = '151016/151016UVVis.csv'
    fluorfolder = '151016/151016fluor/'
    sample1.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dots1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = (406,474)))
    sample2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dots2N2.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = (406,474)))
    sample3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dots1air.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = (406,474)))
    sample4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = (406,474)))
    sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    halfwidths.append(HWHM(loadtxt(uvvisfile,unpack = True, delimiter = ',', skiprows = 1, usecols = (0,1)), (408,420)))
    fr = (474,680)
    #sample1.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dots1N2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = (406,474)))
    midgapemission.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dots1N2.dat', anthracenefile,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = fr))
    mg2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dots2N2.dat', anthracenefile,label = 'N2pH9',color = 'r',fluorescencerange = fr))
    mg3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dots1air.dat', anthracenefile,label = 'N2pH11',color = 'r',fluorescencerange = fr))
    mg4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,label = 'airpH7',color = 'r',fluorescencerange = fr))
    #sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
          
     
     ###################################################
    anthracenefile = '151017/151017fluor/anthracene.dat'
    
    uvvisfile = '151017/151017UVvis.csv'
    fluorfolder = '151017/151017fluor/'
    sample1.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dot1n2.dat', anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = (406,474)))
    sample2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dot2n22.dat', anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = (406,474)))
    sample3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air2.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = (406,474)))
    sample4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r',fluorescencerange = (406,474)))
    sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    halfwidths.append(HWHM(loadtxt(uvvisfile,unpack = True, delimiter = ',', skiprows = 1, usecols = (0,1)), (408,420)))
    
    fr= (474,680)
    midgapemission.append(indivQY(uvvisfile, 'b','g',fluorfolder+'/dot1n2.dat', anthracenefile,day = 0,label = 'N2pH7',color = 'r',fluorescencerange = fr))
    mg2.append(indivQY(uvvisfile, 'c','g',fluorfolder+'/dot2n22.dat', anthracenefile,day = 0,label = 'N2pH9',color = 'r',fluorescencerange = fr))
    mg3.append(indivQY(uvvisfile, 'd','g',fluorfolder+'/dot1air2.dat', anthracenefile,day = 0,label = 'N2pH11',color = 'r',fluorescencerange = fr))
    mg4.append(indivQY(uvvisfile, 'e','g',fluorfolder+'/dot2air.dat', anthracenefile,day = 0,label = 'airpH7',color = 'r',fluorescencerange = fr))
    #sample5.append(indivQY(uvvisfile, 'f','g',fluorfolder+'/dotswithKCl.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r',fluorescencerange = (406,474)))
    
    
    n2samples = mean([array(sample1)/sample1[0],array(sample2)/sample2[0]],axis = 0)
    airsamples = mean([array(sample3)/sample3[0],array(sample4)/sample4[0]],axis = 0)
    
    figure()
    plot(sample1,'ko-',label = 'dot1N2')
    plot(midgapemission,'k--o',label='midgap1N2')
    plot(sample2,'ks-',label = 'dot2N2')
    plot(sample3,'ro-',label = 'dot1air')
    plot(mg3,'r--o',label='mg1air')
    plot(sample4,'rs-',label = 'dot2air')
   # plot(sample5,'bs-',label = 'dots with KCl')
    plot(mg2,'k--s',label='midgap2N2')
    plot(mg4,'r--s',label='midgap2air')
       
    legend()
    ylabel('band edge QY')
    xlabel('day')
    ylim(0,0.0009)
    
    fig4 = figure()
    plot(n2samples,label='N$_2')
    plot(airsamples,label='air')
    def exponentialdecay(x,r,A0,Ainf):return (A0-Ainf)*exp(r*x)+Ainf
    fit = curve_fit(exponentialdecay,arange(0,6),n2samples,[-1,0.7,0.3])
    plot(arange(0,6,0.1),exponentialdecay(arange(0,6,0.1),*fit[0]),'r--')
    print 'n2',fit[0],np.sqrt(np.diag(fit[1]))
    
    fit = curve_fit(exponentialdecay,arange(0,6),airsamples,[-1,0.7,0.3])
    plot(arange(0,6,0.1),exponentialdecay(arange(0,6,0.1),*fit[0]),'k--')
    print 'air',fit,np.sqrt(np.diag(fit[1]))
    legend()
    ylabel('Normalized band edge QY')
    xlabel('day')
    ylim(0,1.1)
    Aax1.legend(list((str(i) for i in range(6))))
    i = 1
#    for l in fig4.axes[0].lines:
#        
#        name = 'stabilityfigure'+str(i)+'.txt'
#        savetxt('/home/chris/Desktop/'+name,transpose([l.get_xdata(),l.get_ydata()]))
#        i+=1
#    
#    for l in ax2.lines:
#        name = 'stabilityfigurefluorescence'+str(i)+'.txt'
#        savetxt('/home/chris/Desktop/'+name,transpose([l.get_xdata(),l.get_ydata()]))
#        i+=1
#    for l in Aax4.lines:
#        xx = l.get_xdata()
#        x1 = argmin(abs(xx-480))
#        
#        yy=l.get_ydata()
#        yy[:]-=yy[x1]
#        x = findpeak(xx,yy,(400,420))
#        B=x[0]
#        d=-0.000000066521*B**3+0.00019557*B**2-0.092352*B+13.29
#        eps =21536*d**2.3
#        print 'conc', x[1]/eps
#        name = 'stabilityfigureabsorbance'+str(i)+'.txt'
#        savetxt('/home/chris/Desktop/'+name,transpose([l.get_xdata(),l.get_ydata()]))
#        i+=1
   
    

    

    print halfwidths
    return (sample1, sample2, sample3, sample4, sample5, )



universalfilelist = [ 'dot1d_N2pH7', 'dot1d_N2pH9', 'dot1d_N2pH11', 
                    'dot1d_airpH7', 'dot1d_airpH9', 'dot1d_airpH11',
                    'dot2d_N2pH7', 'dot2d_N2pH9','dot2d_N2pH11',
                    'dot2d_airpH7', 'dot2d_airpH9','dot2d_airpH11' ]
                    
O2_quench_correction = 1+0.00145*158

def correctionfactorforunfilledcuvettes():
#     print 'The data from september6 was incorrect, since the cuvettes were not completely full, so that the fluorescence from the'
#     print 'samples was too high.  Below are the points collected for Sept7 in the "high" position (so that liquid covered the whole probe area"'
#     print 'and "low" such that the sample was in the same position as Sept6.'
#     print 'the ratio of "hi"to"low" is'
     os.chdir('150907/150907fluor')
    
     
     filelist = ['dot1d1N2pH11h9', 'dot1d1N2pH11l8', 
                    'dot1d1N2pH7hi1','dot1d1N2pH7lo1',
                    'dot1d1N2pH7hi4', 'dot1d1N2pH7lo2',
                    'dot1d1N2pH9hi7', 'dot1d1N2pH9lo6',
                    'dot1d1aipH11h6', 'dot1d1aipH11l5',
                    'dot1d1aipH7h1', 'dot1d1aipH7l1', 
                    'dot1d1aipH9h4', 'dot1d1aipH9l3', 
                    'dot2d1N2pH11h5', 'dot2d1N2pH11l4',
                    'dot2d1N2pH7h1', 'dot2d1N2pH7l10', 
                    'dot2d1N2pH9h1', 'dot2d1N2pH9l1', 
                    'dot2d1aipH11h6', 'dot2d1aipH11l5',
                    'dot2d1aipH7h1', 'dot2d1aipH7l1',
                    'dot2d1aipH9h4', 'dot2d1aipH9l1']
     
     thelist = list()
     for i in range(len(filelist)/2):
         if filelist[2*i][0:11]!=filelist[2*i+1][0:11]:
             print filelist[2*i],filelist[2*i+1]
             continue
         else:
             a = loadtxt(filelist[2*i],delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
             hi = RamanSpectrum(pandas.Series(a[1],a[0]))
             b = loadtxt(filelist[2*i+1],delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
             lo = RamanSpectrum(pandas.Series(b[1],b[0]))
             hi.plot()
             lo.plot()
             hidivlo= hi.calc_area((410,473))/lo.calc_area((410,473))
             thelist.append(hidivlo)
     thelist=array(thelist)
#     print mean(thelist)
#     print 'with a standard deviation of ', std(thelist)
#     print  "The values are returned here"
     ## DATA FROM SEPT6 SHOULD BE MULTIPLIED BY THIS NUMBER TO GET THE CORRECT VALUE
     return numpy.mean(thelist)
#     


def aminoanthracene_fluorescenceyield_determination():
    amino = loadtxt('150929/150929fluor/aminoanthracene.dat',delimiter = '\t', unpack = True, usecols=(0,3),skiprows = 1)
    anthracene = loadtxt('150929/150929fluor/anthracene.dat',delimiter = '\t', unpack = True, usecols=(0,3),skiprows = 1)
    amino[1]-=min(amino[1])
    anthracene[1]-=min(anthracene[1])   
    x = argmin(abs(450-amino[0]))
    Iamino = amino[1,x]
    
    x = argmin(abs(420-anthracene[0]))
    Ianth_420 = anthracene[1,x]
    
    plot(amino[0],amino[1])
    plot(anthracene[0],anthracene[1])

    uv = loadtxt('150929/150929UVVis.csv', skiprows = 1, unpack=True, delimiter = ',',usecols=(0,7,8))
    Absspecanth = RamanSpectrum(pandas.Series(uv[2][::-1],uv[0][::-1]))
    Absspecamino = RamanSpectrum(pandas.Series(uv[1][::-1],uv[0][::-1]))
    
    Absspecamino-=min(Absspecamino)
    Absspecanth-=min(Absspecanth)
    Aanth = Absspecanth[350]
    Aamino = Absspecamino[350]
    
    O2quench_correction = (1+0.00145*158)
    
    QY = 0.27/O2quench_correction*Iamino*(1-10**(-Aanth))/(Ianth_420*76.9811)/(1-10**(-Aamino))
    print 'dQ/dlambda for aminoanthracene at 450 is',  QY
    figure()
    ax=subplot(111)
    r = indivQY('150929/150929UVVis.csv','h','i', '150929/150929fluor/aminoanthracene.dat','150929/150929fluor/anthracene.dat',fluorescencerange = (450,455), UVVisplot=ax, fluorplot = ax )
      
    return r
    

    

def trial1_Sept6through13():
    """Photoluminescence of CdS QDs PPA capped. """
    sample1=list()
    sample2=list()
    sample3=list()
    sample4=list()
    sample5=list()
    sample6=list()
    sample7=list()
    sample8=list()
    sample9=list()
    sample10=list()
    sample11=list()
    sample12=list()
    
    figure()
    Aax1 = subplot(231)
    Aax2 = subplot(232)
    Aax3 = subplot(233)
    Aax4 = subplot(234)
    Aax5 = subplot(235)
    Aax6 = subplot(236)
    
    figure()
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(234)
    ax5 = subplot(235)
    ax6 = subplot(236)
    ##day0
  
    uvvisfile = '150906/150906.csv'  #, delimiter = ',',usecols=[0,4,5,6,9,10,11,14,15,16,19,20,21,23]
    anthracenefile = '150906/150906fluor/aminoanthracene.dat'
    fluorfolder = '150906/150906fluor'
    sample1.append(0.907*indivQY(uvvisfile, 'e','i',fluorfolder+'/dot1_O2_pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0.907*indivQY(uvvisfile, 'f','i',fluorfolder +'/dot1_O2_pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(0.907*indivQY(uvvisfile, 'g','i',fluorfolder +'/dot1_O2_pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(0.907*indivQY(uvvisfile, 'j','i',fluorfolder +'/dot1_N2_pH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(0.907*indivQY(uvvisfile, 'k','i',fluorfolder +'/dot1_N2_pH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(0.907*indivQY(uvvisfile, 'l','i',fluorfolder +'/dot1_N2_pH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(0.907*indivQY(uvvisfile, 'o','i',fluorfolder+'/dot2_O2_pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(0.907*indivQY(uvvisfile, 'p','i',fluorfolder +'/dot2_O2_pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(0.907*indivQY(uvvisfile, 'q','i',fluorfolder +'/dot2_O2_pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(0.907*indivQY(uvvisfile, 't','i',fluorfolder +'/dot2_N2_pH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(0.907*indivQY(uvvisfile, 'u','i',fluorfolder +'/dot2_N2_pH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(0.907*indivQY(uvvisfile, 'v','i',fluorfolder +'/dot2_N2_pH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    
    ##day1
    uvvisfile = '150907/150907UVVis.csv'  #, delimiter = ',',usecols=[0,4,5,6,9,10,11,14,15,16,19,20,21,23]
    anthracenefile = '150907/150907fluor/aminoanthracene.dat'
    fluorfolder = '150907/150907fluor/'
    sample1.append(indivQY(uvvisfile, 'b','p',fluorfolder+'/dot1d1N2pH7hi1',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
   # sample2.append(indivQY(uvvisfile, 'c','p',fluorfolder + '/dot1d1N2pH9hi7',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','p',fluorfolder +'/dot1d1N2pH11h9', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','p',fluorfolder +'/dot1d1aipH7h1', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','p',fluorfolder +'/dot1d1aipH9h4', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','p',fluorfolder +'/dot1d1aipH11h6',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','p',fluorfolder+'/dot2d1N2pH7h1',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','p',fluorfolder + '/dot2d1N2pH9h1',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','p',fluorfolder +'/dot2d1N2pH11h5', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','p',fluorfolder +'/dot2d1aipH7h1', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','p',fluorfolder +'/dot2d1aipH9h4', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','p',fluorfolder +'/dot2d1aipH11h6',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
  
    ##day2

     
    uvvisfile = '150908/PPAstability.csv'
    anthracenefile = '150908/150908fluor/aminoanthracene.dat'
    fluorfolder = '150908/150908fluor/'
    sample1.append(indivQY(uvvisfile, 'b','o',fluorfolder+'/dot1d2N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0)#ivQY(uvvisfile, 'c','o',fluorfolder + '/dot1d2N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','o',fluorfolder +'/dot1d2N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','o',fluorfolder +'/dot1d2airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','o',fluorfolder +'/dot1d2airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','o',fluorfolder +'/dot1d2airpH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','o',fluorfolder+'/dot2d2N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','o',fluorfolder + '/dot2d2N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','o',fluorfolder +'/dot2d2N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','o',fluorfolder +'/dot2d2airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','o',fluorfolder +'/dot2d2airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','o',fluorfolder +'/dot2d2airpH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
   
    ##day3
    uvvisfile = '150909/150909UVVis.csv'
    anthracenefile = '150909/150909fluor/aminoanthracene.dat'
    fluorfolder = '150909/150909fluor/'
    sample1.append(indivQY(uvvisfile, 'b','o',fluorfolder+'/dot1d3N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
   # sample2.append(0)#ivQY(uvvisfile, 'c','o',fluorfolder + '/dot1d2N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','o',fluorfolder +'/dot1d3N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','o',fluorfolder +'/dot1d3airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','o',fluorfolder +'/dot1d3airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','o',fluorfolder +'/dot1d3airpH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','o',fluorfolder+'/dot2d3N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','o',fluorfolder + '/dot2d3N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','o',fluorfolder + '/dot2d3N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','o',fluorfolder +'/dot2d3airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','o',fluorfolder +'/dot2d3airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','o',fluorfolder +'/dot2d3airpH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
#    ##day4
 
    uvvisfile ='150910/150910UVVis.csv'
    anthracenefile = '150910/150910fluor/aminoanthracene.dat'
    fluorfolder = '150910/150910fluor/'
    sample1.append(indivQY(uvvisfile, 'b','n',fluorfolder+'/dot1d4N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
   # sample2.append(0)#ivQY(uvvisfile, 'c','n',fluorfolder + '/dot1d2N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','n',fluorfolder +'/dot1d4N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','n',fluorfolder +'/dot1d4airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','n',fluorfolder +'/dot1d4airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','n',fluorfolder +'/dot1d4airpH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','n',fluorfolder+'/dot2d4N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','n',fluorfolder +'/dot2d4N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','n',fluorfolder +'/dot2d4N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','n',fluorfolder +'/dot2d4airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','n',fluorfolder +'/dot2d4airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','n',fluorfolder +'/dot2d4airpH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
#    ##day5
#    day=5
    uvvisfile ='150911/150911UVVis.csv'
    anthracenefile = '150911/150911fluor/aminoanthracene.dat'
    fluorfolder = '150911/150911fluor/'
    sample1.append(indivQY(uvvisfile, 'b','o',fluorfolder+'/dot1d5N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0)#ivQY(uvvisfile, 'c','o',fluorfolder + '/dot1d2N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','o',fluorfolder +'/dot1d5N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','o',fluorfolder +'/dot1d5airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','o',fluorfolder +'/dot1d5airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','o',fluorfolder +'/dot1d5airpH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','o',fluorfolder+'/dot2d5N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','o',fluorfolder + '/dot2d5N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','o',fluorfolder +'/dot2d5N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','o',fluorfolder +'dot2d5airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','o',fluorfolder +'dot2d5airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','o',fluorfolder +'dot2d5airpH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
#     
#      ##day7
    day=7
    uvvisfile ='150913/150913UVVis.csv'
    anthracenefile = '150913/150913fluor/aminoanthracene.dat'
    fluorfolder = '150913/150913fluor/'
    sample1.append(0)#indivQY(uvvisfile, 'b','n',fluorfolder+'dot1d7N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0)#ivQY(uvvisfile, 'c','n',fluorfolder + 'dot1d2N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','n',fluorfolder +'dot1d7N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','n',fluorfolder +'dot1d7airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','n',fluorfolder +'dot1d7airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','n',fluorfolder +'dot1d7airpH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','n',fluorfolder+'dot2d7N2pH7.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','n',fluorfolder + 'dot2d7N2pH9.dat',  anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','n',fluorfolder +'dot2d7N2pH11.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','n',fluorfolder +'dot2d7airpH7.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','n',fluorfolder +'dot2d7airpH9.dat', anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','n',fluorfolder +'dot2d7airpH11.dat',anthracenefile, fluorescencerange = (408,480),UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
    
    figure()
    ax1 = subplot(131)
    ax2 = subplot(132)
    ax3=subplot(133)
    days = [0,1,2,3,4,5,7]
    ax1.plot(days,sample1,'rs-',label = 'N2pH7')
   # plot(days,sample2,'rs-',label = 'N2pH9')
    ax3.plot(days,sample3,'rs-',label = 'N2pH11')
    ax1.plot(days,sample4,'ks-',label = 'airpH7')
    ax2.plot(days,sample5,'ks-',label = 'airpH9')
    ax3.plot(days,sample6,'ks-',label = 'airpH11')
    ax1.plot(days,sample7,'ro-',label = 'N2pH7')
    ax2.plot(days,sample8,'ro-',label = 'N2pH9')
    ax3.plot(days,sample9,'ro-',label = 'N2pH11')
    ax1.plot(days,sample10,'ko-',label = 'airpH7')
    ax2.plot(days,sample11,'ko-',label = 'airpH9')
    ax3.plot(days,sample12,'ko-',label = 'airpH11')
    ax1.set_ylim(0,0.00075)
    ax2.set_ylim(0,0.00075)
    ax3.set_ylim(0,0.00075)
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax1.set_ylabel('band edge QY')
    ax1.set_xlabel('day')
    ax2.set_xlabel('day')
    ax3.set_xlabel('day')
    return (sample1, sample2, sample3, sample4, sample5, sample6,sample7, sample8,sample9,sample10,sample11,sample12)   