# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 15:13:15 2016

@author: chris
"""
from ramanTools.OrcaTools import NWChemDOS 
import numpy as np
from ramanTools.RamanSpectrum import *
from ramanTools.OrcaTools import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import numpy as np

import pickle
from makeacluster import pointsonsphere

import load_cube
from math import factorial as fac
from UVVistools import *

def calculateFrohlichcoupling(n,m):
    """Calculate the Frohlich coupling from resonance Raman spectra.  See alivisatos paper on resonance Raman of CdSe"""
    if material == 'CdS':
        m_e = 0.18
        m_h = 0.51
        a0 = 5.82E-10   ## Lattice constant (m)
        eps_bulk = 5.3
        eps_0 = 8.7
        bohr = 21.6 # bohr exciton radius (m)
        wLO = 2*pi*1240/(0.01/305)/h  ### s-1
        delta= 3.07
    elif material == 'CdSe':
        m_e = 0.13
        m_h = 0.45
        a0 = 6.05E-10 ## Lattice constant (m)
        eps_bulk = 6.1
        eps_0 = 9.3
        bohr = 32.3  # bohr exciton radius (m)
        wLO = 2*pi*1240/(0.01/210)/h  ### s-1
        delta =  2.93
        
    e=1.602E-19  # C
       
    h = 6.626E-36  #Js
    wLO = 2*pi*1240/(0.01/205)/h
    print wLO
    eps_bulk = 10.16
    eps_0 = 1
    bohr=1
    w = (3*pi**2)**(1/3)*(a0/bohr)
    
    def L(i,j,x):
        return 1
    
    x = linspace(0,w,1000)
    y = x**4*(2+x**2)**2*(1+x**2)**-4
    sumy = sum(y)*w/1000
    
    delta = 1.97*e**2/(a0*h*wLO)*(1/eps_bulk-1/eps_0)*(1/w)  * sumy
    
    densitymatrix = np.sqrt(fac(n)/fac(m))*exp(-0.5*delta**2)*delta**(n-m)*L(m,n-m,delta**2)
    return densitymatrix



def CdPPAReferences():
    """non-resonant Raman spectra of Cd-PPA complexes from January 19, 2016"""
    
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160119/160119_01.txt')  ## pH 9
    b= RamanSpectrum('/home/chris/Dropbox/DataWeiss/160119/160119_02.txt')  ## pH 11
    c= RamanSpectrum('/home/chris/Dropbox/DataWeiss/160119/160119_03.txt')  ## pH 12
    d= RamanSpectrum('/home/chris/Dropbox/DataWeiss/160119/160119_04.txt')  ## pH 13
    e= RamanSpectrum('/home/chris/Dropbox/DataWeiss/160113/160113_04.txt')  ## pH 5
    f = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160113/160113_05.txt')  ## pH 6
    g = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160113/160113_06.txt')  ## pH 8
    h= RamanSpectrum('/home/chris/Dropbox/DataWeiss/160113/160113_03.txt') #CdPPA\

    offset=0
    for x in [e,f,g,a,b,c,d,h]:
        x.smooth()
        x.autobaseline((200,1700), order = 4)
        x.autobaseline((1700,4000),order = 4, join='start')
        x[:]/=max(x[200:1700])
        x[:]+=offset
        offset+=1
        x.plot()    
    xlim(200,3800)
    annotate('pH5',(2000,0.2),fontsize = 20)
    annotate('pH6',(2000,1.2),fontsize = 20)
    annotate('pH8',(2000,2.2),fontsize = 20)
    annotate('pH9',(2000,3.2),fontsize = 20)
    annotate('pH11',(2000,4.2),fontsize = 20)
    annotate('pH12',(2000,5.2),fontsize = 20)
    annotate('pH13',(2000,6.2),fontsize = 20)
    annotate('Cd(OH)$_2$',(2000,7.2),fontsize = 20)

    
    
    
    return 0
    
def Jan19():  ### PPA on QDs crashed from hexanes.
    """Raman spectra of CdS dots at different steps in PPA exchange January 19"""

    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160119/160119_05.txt')
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160119/160119_06.txt')
    c = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160119/160119_08.txt')

    offset=0
    for x in [a,b,c]:
        x.smooth()
        x.autobaseline((200,1700), order = 4)
        #x.autobaseline((1700,4000),order = 4, join='start')
        x[:]/=max(x[200:1700])
        x[:]+=offset
        offset+=1
        x.plot()
    xlim(200,3800)
    return 0
    

def Jan11():
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160111/160111_03.txt')
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160111/160111_04.txt')
    c = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160111/160111_02.txt')
    
    offset=0
    for x in [a,b,c]:
        print x.name
        x.smooth()
        x.autobaseline((200,1700), order = 4)
        #x.autobaseline((1700,4000),order = 4, join='start')
        x[:]/=max(x[200:1700])
        x[:]+=offset
        offset+=1
        x.plot()
    return 0

def Jan6():  ## CdPPA dots through the exchange
    """CdPPA dots through the exchange"""
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_03.txt')
    b= RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_05.txt')
    c= RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_06.txt')
    offset=0
    for x in [a,b,c]:
        print x.name
        x.smooth()
        x.autobaseline((200,1800), order = 4)
        #x.autobaseline((1700,4000),order = 4, join='start')
        x[:]/=max(x[200:1800])
        x[:]+=offset
        offset+=1
        x.plot()
    return 0
    
    

    
def CdSeOleateSynthJan29():
    """View the UV vis spectrum of oleate capped CdSe quantum dots synthesiszed on January 29 2016"""
    
    a = loadtxt('/home/chris/Dropbox/DataWeiss/160129/cdse oleate synth jan29.csv', delimiter = ',', unpack = True, skiprows = 1)
    a[1]-=a[1][0]
    b = findpeak(a[0],a[1],(520,545))
    print b
    plot(*a)
    print a[0][argmin(abs(b[1]/2-a[1]))]
    print 'fwhm', a[0][argmin(abs(b[1]/2-a[1]))]-b[0]
    
    concentration = CdSconc(b)
    return concentration



    
def Feb1():  ## 
    """Resonance Raman of CdSe dots with PPA in water.  February 1"""
    clf()
    a473 =  RamanSpectrum('/home/chris/Dropbox/DataWeiss/160201/160201_07.txt')
    a633 =  RamanSpectrum('/home/chris/Dropbox/DataWeiss/160201/160201_08.txt')
    a633=SPIDcorrect633(a633)
    a473.autobaseline((120,700),order = 3)
    a633.autobaseline((120,700),order = 3)
    a473.plot()
    a633.plot()
    legend(['473','633'])
    return 0
    
def Feb3():  
    """UVVis spectra of PPA capped CdS dots in varying concentrations and pH of PPA/KOH buffer"""
   
   
    pHs = array([6.2,    6.5,    7,    7.5,    8,   8.7,
                6.7,    6.9,    7.3,    7.7,    8.1,    8.8,
                7.3,    7.5,    7.6,    8.0,    8.3,    8.8,
                6,      6.5,    7.0,    7.5,    8,    8.6,
                6.7,    7.0,    7.3,    7.7,    8.1,    8.7,
                7.6,   7.6,    7.4,    7.7,    8.1,    8.6,
                6,    6.5,    7,    7.5,    8,    8.6,
                6.8,    6.9,    7.2,    7.6,    8.0,    8.6,
                7.5,    7.4,    7.5,    7.8,    8.1,    8.6]).reshape(3,3,6)
    
    a = loadtxt('/home/chris/Dropbox/DataWeiss/160203/160203_PPAdotsinbuffers.csv', delimiter =',', unpack = True, skiprows = 2)
    
    a= a[[6,]+range(9,116,2)]
    print a.shape
    a[1:]-=a[1,0]
    peaks = array([])
    absorbances = array([])
    for i in a[1:]:
        p = findpeak(a[0],i,(406,415))
        peaks = append(peaks,p[0])
        absorbances = append(absorbances, p[1])
    peaks = peaks.reshape(3,3,-1)
    absorbances = absorbances.reshape(3,3,-1)

    fig1 = figure()
    
    CdS1 = 0
    CdS2 = 1
    CdS3 =2
    mM18=0  
    mM10 = 1 
    mM2= 2
    ##### samples CdS1 with concentrations, 2mM, 10mM, 18mM PPA buffer.  pH vs wavelength
    plot(pHs[CdS1, mM18,:],peaks[CdS1, mM18,:],'rs-')  ### 18 mM 
    plot(pHs[CdS1,1,:],peaks[CdS1,1,:],'bs-')  ### 10 mM 
    plot(pHs[CdS1,2,:],peaks[CdS1,2,:],'ks-')  ## 2 mM
    ##### samples  CdS2, 2mM, 10mM, 18mM PPA buffer.  pH vs wavelength
    plot(pHs[1,0,:],peaks[1,0,:],'ro-')
    plot(pHs[1,1,:],peaks[1,1,:],'bo-')
    plot(pHs[1,2,:],peaks[1,mM2,:],'ko-')
    
    ##### samples  CdS3, 2mM, 10mM, 18mM PPA buffer.  pH vs wavelength
    plot(pHs[2,0,:],peaks[2,0,:],'rx-')
    plot(pHs[2,1,:],peaks[2,1,:],'bx-')
    plot(pHs[2,2,:],peaks[2,mM2,:],'kx-')    
    
    
    
    print '----18 mM PPA slope lambdamaxvspH-----'
    a= numpy.polyfit(pHs[CdS1, mM18,:],peaks[CdS1, mM18,:],1) 
    b= numpy.polyfit(pHs[CdS2, mM18,:],peaks[CdS2, mM18,:],1)
    c= numpy.polyfit(pHs[CdS3, mM18,:],peaks[CdS3, mM18,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)
    print '----10 mM PPA slope lambdamaxvspH-----'
    a= numpy.polyfit(pHs[CdS1, mM10,:],peaks[CdS1, mM10,:],1) 
    b= numpy.polyfit(pHs[CdS2, mM10,:],peaks[CdS2, mM10,:],1)
    c=  numpy.polyfit(pHs[CdS3, mM10,:],peaks[CdS3, mM10,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)
    print '----2 mM PPA slope lambdamaxvspH-----'
    a =  numpy.polyfit(pHs[CdS1, mM2,:],peaks[CdS1, mM2,:],1) 
    b =  numpy.polyfit(pHs[CdS2, mM2,:],peaks[CdS2, mM2,:],1)
    c =  numpy.polyfit(pHs[CdS3, mM2,:],peaks[CdS3, mM2,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)  


    figure()
    peakseV = 1240000/peaks-1240000/408  ### 408 is the wavelength of the oleate capped dots
     ##### samples CdS1 with concentrations, 2mM, 10mM, 18mM PPA buffer.  pH vs wavelength
#    plot(pHs[CdS1, mM18,:],peakseV[CdS1, mM18,:],'rs-')  ### 18 mM 
#    plot(pHs[CdS1,1,:],peakseV[CdS1,1,:],'bs-')  ### 10 mM 
#    plot(pHs[CdS1,2,:],peakseV[CdS1,2,:],'ks-')  ## 2 mM
    
    ##### samples  CdS2, 2mM, 10mM, 18mM PPA buffer.  pH vs wavelength
#    plot(pHs[1,0,:],peakseV[1,0,:],'ro-')
#    plot(pHs[1,1,:],peakseV[1,1,:],'bo-')
#    plot(pHs[1,2,:],peakseV[1,mM2,:],'ko-')
    errorbar(pHs[1,0,:],np.mean(peakseV[:,mM2,:],axis=0),yerr=np.std(peakseV[:,mM2,:]))
    plot(pHs[1,0,:],np.mean(peakseV[:,mM10,:],axis=0))
    plot(pHs[1,0,:],np.mean(peakseV[:,mM18,:],axis=0))
    ylabel('$\Delta\lambda$ first exciton peak from oleate peak (meV)', fontsize=20)
    xlabel('pH', fontsize=20)
    legend(['2 mM PPA', '10 mM PPA', '18 mM PPA'])
    ##### samples  CdS3, 2mM, 10mM, 18mM PPA buffer.  pH vs wavelength
#    plot(pHs[2,0,:],peakseV[2,0,:],'rx-')
#    plot(pHs[2,1,:],peakseV[2,1,:],'bx-')
#    plot(pHs[2,2,:],peakseV[2,mM2,:],'kx-')    
    
    print '----18 mM PPA slope lambdamaxvspH-----'
    a= numpy.polyfit(pHs[CdS1, mM18,:],peakseV[CdS1, mM18,:],1) 
    b= numpy.polyfit(pHs[CdS2, mM18,:],peakseV[CdS2, mM18,:],1)
    c= numpy.polyfit(pHs[CdS3, mM18,:],peakseV[CdS3, mM18,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)
    print '----10 mM PPA slope lambdamaxvspH-----'
    a= numpy.polyfit(pHs[CdS1, mM10,:],peakseV[CdS1, mM10,:],1) 
    b= numpy.polyfit(pHs[CdS2, mM10,:],peakseV[CdS2, mM10,:],1)
    c=  numpy.polyfit(pHs[CdS3, mM10,:],peakseV[CdS3, mM10,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)
    print '----2 mM PPA slope lambdamaxvspH-----'
    a =  numpy.polyfit(pHs[CdS1, mM2,:],peakseV[CdS1, mM2,:],1) 
    b =  numpy.polyfit(pHs[CdS2, mM2,:],peakseV[CdS2, mM2,:],1)
    c =  numpy.polyfit(pHs[CdS3, mM2,:],peakseV[CdS3, mM2,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)  
    
    figure()
    suptitle('absorbance')
    
     ##### samples CdS1 with concentrations, 2mM, 10mM, 18mM PPA buffer.  pH vs absorbance
    plot(pHs[CdS1, mM18,:],absorbances[CdS1, mM18,:],'rs-')  ### 18 mM 
    plot(pHs[CdS1,1,:],absorbances[CdS1,1,:],'bs-')  ### 10 mM 
    plot(pHs[CdS1,2,:],absorbances[CdS1,2,:],'ks-')  ## 2 mM
    ##### samples  CdS2, 2mM, 10mM, 18mM PPA buffer.  pH vs absorbance
    plot(pHs[1,0,:],absorbances[1,0,:],'ro-')
    plot(pHs[1,1,:],absorbances[1,1,:],'bo-')
    plot(pHs[1,2,:],absorbances[1,mM2,:],'ko-')
    
    ##### samples  CdS3, 2mM, 10mM, 18mM PPA buffer.  pH vs absorbance
    plot(pHs[2,0,:],absorbances[2,0,:],'rx-')
    plot(pHs[2,1,:],absorbances[2,1,:],'bx-')
    plot(pHs[2,2,:],absorbances[2,mM2,:],'kx-')    
    
    print '----18 mM PPA slope absorbancevspH-----'
    a= numpy.polyfit(pHs[CdS1, mM18,:],absorbances[CdS1, mM18,:],1) 
    b= numpy.polyfit(pHs[CdS2, mM18,:],absorbances[CdS2, mM18,:],1)
    c= numpy.polyfit(pHs[CdS3, mM18,:],absorbances[CdS3, mM18,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)
    print '----10 mM PPA slope absorbancevspH-----'
    a= numpy.polyfit(pHs[CdS1, mM10,:],absorbances[CdS1, mM10,:],1) 
    b= numpy.polyfit(pHs[CdS2, mM10,:],absorbances[CdS2, mM10,:],1)
    c=  numpy.polyfit(pHs[CdS3, mM10,:],absorbances[CdS3, mM10,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)
    print '----2 mM PPA slope absorbancevspH-----'
    a =  numpy.polyfit(pHs[CdS1, mM2,:],absorbances[CdS1, mM2,:],1) 
    b =  numpy.polyfit(pHs[CdS2, mM2,:],absorbances[CdS2, mM2,:],1)
    c =  numpy.polyfit(pHs[CdS3, mM2,:],absorbances[CdS3, mM2,:],1)
    datas = array([a[0],b[0],c[0]])
    print 'average slope', np.mean(datas), 'stdev:', np.std(datas)  
    
    return 0
    
def Feb4():  
    """UVVis spectra of PPA capped CdS dots with partial displacement by MPA. Experiment Feb 4"""
    
    pHs = array([6.3,    7,    8,   8.6]*18).reshape(3,6,4)
    pHs[0,0,0] = 6.0
    a = loadtxt('/home/chris/Dropbox/DataWeiss/160204/PPAcappeddotswithMPAdisplacement.csv', delimiter =',', unpack = True, skiprows = 2)
    
    a= a[[0,]+range(1,144,2)]
    print a.shape
    a[1:]-=a[1,0]
    peaks = array([])
    absorbances = array([])
    for i in a[1:]:
        
        p = findpeak(a[0],i,(405,413))
        peaks = append(peaks,p[0])
        absorbances = append(absorbances, p[1])
    
    peaks = peaks.reshape(6,3,4)
    
    peaks = np.transpose(peaks,axes = (1,0,2))
    absorbances = absorbances.reshape(6,3,4)
    absorbances = np.transpose(absorbances, axes =(1,0,2) )
    
    fig1 = figure()
    
    CdS1 = 0
    CdS2 = 1
    CdS3 =2
    MPA0 = 0
    MPA50=1
    MPA100 = 2
    MPA150 = 3
    MPA200 =4
    MPA250 = 5
    
  
    #### samples CdS1 with concentrations, 0,50,100,150,200,250 equivalents MPA. pH vs wavelength
    plot(pHs[CdS1, MPA0,:],peaks[CdS1, MPA0,:],'rs-')  ### 0 equivalents pH vs wavelength
    plot(pHs[CdS1,MPA50,:],peaks[CdS1,MPA50,:],'bs-')  ### 50 eq
    plot(pHs[CdS1,MPA100,:],peaks[CdS1,MPA100,:],'gs-')  ## 100 eq
    plot(pHs[CdS1,MPA150,:],peaks[CdS1,MPA150,:],'ys-')  ## 150 eq
    plot(pHs[CdS1,MPA200,:],peaks[CdS1,MPA200,:],'cs-')  ## 200 eq
    plot(pHs[CdS1,MPA250,:],peaks[CdS1,MPA250,:],'ks-')  ## 250eq 
    ##### samples  CdS2, 2mM, 10mM, 18mM PPA buffer.  pH vs wavelength
    plot(pHs[CdS2, MPA0,:],peaks[CdS2, MPA0,:],'ro-')  ### 0 equivalents pH vs wavelength
    plot(pHs[CdS2,MPA50,:],peaks[CdS2,MPA50,:],'bo-')  ### 50 eq
    plot(pHs[CdS2,MPA100,:],peaks[CdS2,MPA100,:],'go-')  ## 100 eq
    plot(pHs[CdS2,MPA150,:],peaks[CdS2,MPA150,:],'yo-')  ## 150 eq
    plot(pHs[CdS2,MPA200,:],peaks[CdS2,MPA200,:],'co-')  ## 200 eq
    plot(pHs[CdS2,MPA250,:],peaks[CdS2,MPA250,:],'ko-')   ## 250eq 
    
    plot(pHs[CdS3, MPA0,:],peaks[CdS3, MPA0,:],'rx-')  ### 0 equivalents pH vs wavelength
    plot(pHs[CdS3,MPA50,:],peaks[CdS3,MPA50,:],'bx-')  ### 50 eq
    plot(pHs[CdS3,MPA100,:],peaks[CdS3,MPA100,:],'gx-')  ## 100 eq
    plot(pHs[CdS3,MPA150,:],peaks[CdS3,MPA150,:],'yx-')  ## 150 eq
    plot(pHs[CdS3,MPA200,:],peaks[CdS3,MPA200,:],'cx-')  ## 200eq
    plot(pHs[CdS3,MPA250,:],peaks[CdS3,MPA250,:],'kx-')  ## 250 eq 

    figure()
    title('averages')
    peaks = 1240000/peaks
    errs = array([])
    slopes = array([])
    for i in [MPA0, MPA50, MPA100, MPA150,  MPA200, MPA250]:
        x= [6.3,    7,    8,   8.6]
        y = np.mean(peaks[:, i,:],axis = 0 )
        if i == MPA150 or i ==MPA0:
            yerr = np.std(peaks[:, i,:],axis = 0 )
        else: yerr = 0
        
        errorbar(x, y, yerr=yerr,marker = 's')
        print '----',i*50,'MPA slope lambdamaxvspH-----'
        a= numpy.polyfit(pHs[CdS1, i,:],peaks[CdS1, i,:],1) 
        b= numpy.polyfit(pHs[CdS2, i,:],peaks[CdS2, i,:],1)
        c= numpy.polyfit(pHs[CdS3, i,:],peaks[CdS3, i,:],1)
        print a[0],b[0],c[0]
        datas = array([a[0],b[0],c[0]])
        print 'average slope', np.mean(datas), 'stdev:', np.std(datas)
        slopes = append(slopes, np.mean(datas))
        errs = append(errs, np.std(datas))
    legend(['0eq MPA','50eq', '100eq','150eq','200eq','250eq'])
    print slopes
    figure()
    
    errorbar(arange(0,300,50), slopes, yerr=errs, marker = 's')
    
    return 0    