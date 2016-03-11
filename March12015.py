# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 17:53:28 2016

@author: chris
"""
import weissdatavariables
from ramanTools.RamanSpectrum import *
import UVVistools
import numpy as np
import ramanTools.RamanTools
import Tkinter
def March1():
    fluorfilelist = ['160301/160301fluor/160301_01.txt',  ### PPA capped
                     '160301/160301fluor/160301_02.txt', ### oleate
                     
                    ]
                     
    filenames = array(['1pH5', '2pH5', '3pH5', '4pH5',
                       '1pH7', '2pH7', '3pH7', '4pH7',
                       '1pH8', '2pH8', '3pH8', '4pH8',
                       '1pH9', '2pH9', '3pH5', '4pH9',
                       '1pH11', '2pH11', '3pH11', '4pH11',])
                     
    uvvis = loadtxt('160301/160301.csv', delimiter = ',', unpack = True,skiprows = 1,usecols=(0,5,6,8))
    uvvis[1:]-=transpose([uvvis[1:,0]])
    
    x448nm = np.where(uvvis[0]==448)[0][0]
    x473nm = np.where(uvvis[0]==473)[0][0]
    print x473nm
    uvvis[3]-=uvvis[3,x448nm]
    
    for i in uvvis[1:]: plot(uvvis[0],i)
    figure()
    absorbances_dots = uvvis[1:3,x473nm]
    plot(absorbances_dots)
    figure()
    absorbance_anth = uvvis[-1,x473nm]
    
    fluorescencedots = array([])
    nliq = 1.333
    nE = 1.359
    
    
    s =  RamanSpectrum('160301/160301fluor/160301_03.txt')
    indy = 10**7/(10**7/473-s.index)
    rhodB = RamanSpectrum(pandas.Series(s.values, indy))
    rhodBarea = rhodB.calc_area((500,700))
    rhodB*=0.65*(1-10**(-absorbance_anth))*nE**2/nE**2 /rhodBarea/(1-10**(-absorbance_anth)) 
    print rhodB.calc_area((500,700))
    #rhodB.plot()
    
    
    
    ## PPA capped dots
    s = RamanSpectrum(fluorfilelist[0])
    indy = 10**7/(10**7/473-s.index)
    s = RamanSpectrum(pandas.Series(s.values, indy))
    s*=0.65*(1-10**(-absorbance_anth))*nliq**2/nE**2 /rhodBarea/(1-10**(-absorbances_dots[0])) 
    s.plot()
    print 'PPA capped CdSe dots QY:', s.calc_area((490,650))
    
    
    #### oleate dots
    s = RamanSpectrum(fluorfilelist[1])
    indy = 10**7/(10**7/473-s.index)
    s = RamanSpectrum(pandas.Series(s.values, indy))
    s*=0.65*(1-10**(-absorbance_anth))*1.375**2/nE**2/rhodBarea/(1-10**(-absorbances_dots[1])) 
    s.plot()
    print 'oleate capped CdSe dots QY:', s.calc_area((490,650))

    
    return 0
    
def March3():
    """Titration of CdS dots PPA capped.  Using PPA as acid"""
    a = loadtxt('160303/160303_titrationpart-allparts.csv', delimiter = ',',skiprows = 2, usecols =(0,)+tuple(range(1,78,2)),unpack = True)
    b = loadtxt('160303/TitrationAdditionsAndpHs.csv',delimiter=',',skiprows = 1, unpack = True)
    pHs = b[5]
    volumes = b[9]
    peakxs = array([])
    peakys = array([])
    halfwidths  = array([])
    a[1:]*=transpose([volumes/volumes[0]])
    a[1:]-=transpose([a[1:,0]])
    figure()
    for i in a[1:]:
#        plot(a[0],i)
#        (x,y) = UVVistools.findpeak(a[0],i,(410,420))
#        peakxs = np.append(peakxs,x)
#        peakys = np.append(peakys,y-i[45])
#        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(410,420)))
        
        (x,y) = UVVistools.findpeak(a[0],i,(401,420))
        peakxs = np.append(peakxs,x)
        z  = argmin(abs(x+50-a[0]))
        i-=i[z]
        peakys = np.append(peakys,y-i[z])
        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(401,420),_plot=False))
        plot(a[0],i)
   
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    ax1.plot(pHs)
    
    ax2.plot(peakxs,'o-')
    ax3.plot(peakys,'o-')
    ax4.plot(halfwidths,'o-')
    
    ax1.set_ylabel('pH')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
    
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    ax2.plot(pHs,peakxs,'o-')
    ax3.plot(pHs,peakys,'o-')
    ax4.plot(pHs,halfwidths,'o-')
    

    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
  
    figure()
    for i in a[1:16]:
        plot(a[0],i)
    return 0    
    
    
        
def March4():
    """Titration of CdS dots PPA capped.  Using PPA as acid.  Second attempt"""
    a = loadtxt('160304/160304_PPACdSTitrationAll.csv', delimiter = ',',skiprows = 2, usecols =(0,)+tuple(range(3,106,2)),unpack = True)
    b = loadtxt('160304/160304_TitrationAdditionsAndpHs.csv',delimiter=',',skiprows = 1, unpack = True)
    pHs = b[5]
    volumes = b[9]
    peakxs = array([])
    peakys = array([])
    halfwidths  = array([])
    a[1:]*=transpose([volumes/volumes[0]])
    a[1:]-=transpose([a[1:,0]])
    
    figure()
    for i in a[1:]:
#        plot(a[0],i)
#        (x,y) = UVVistools.findpeak(a[0],i,(401,418))
#        peakxs = np.append(peakxs,x)
#        peakys = np.append(peakys,y-i[345])
#        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(401,418),_plot=True))
        
        (x,y) = UVVistools.findpeak(a[0],i,(401,418))
        peakxs = np.append(peakxs,x)
        z  = argmin(abs(x+50-a[0]))
        i-=i[z]
        peakys = np.append(peakys,y-i[z])
        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(401,418),_plot=False))
        plot(a[0],i)
   
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    ax1.plot(pHs)
    
    ax2.plot(peakxs,'o-')
    ax3.plot(peakys,'o-')
    ax4.plot(halfwidths,'o-')
    
    ax1.set_ylabel('pH')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
    
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    
    ax1.plot(pHs,peakys*halfwidths)
    ax2.plot(pHs,peakxs,'o-')
    ax3.plot(pHs,peakys,'o-')
    ax4.plot(pHs,halfwidths,'o-')
    
    ax1.set_ylabel('first exicton areas (relative)')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
  
 
    figure()
    for i in a[29:48:2]:
        d =  UVVistools.findpeak(a[0],i,(401,418))[0]
        plot(a[0]-d,i)
    return 0    
    
def March5():
    """Stability of  CdS dots PPA capped. """
    a = loadtxt('160305/160305_stability measurement.csv', delimiter = ',',skiprows = 2, usecols =(0,)+tuple(range(75,110,2)),unpack = True)
    b = loadtxt('160305/160305_TitrationAdditionsAndpHs.csv',delimiter=',',skiprows = 1, unpack = True)
    pHs = b[5]
    volumes = b[9]
    volumes-=0.0002
    peakxs = array([])
    peakys = array([])
    halfwidths  = array([])
    a[1:]*=transpose([volumes/volumes[0]])
    a[1:]-=transpose([a[1:,0]])
    
    figure()
    for i in a[1:]:
        
        (x,y) = UVVistools.findpeak(a[0],i,(401,418))
        peakxs = np.append(peakxs,x)
        z  = argmin(abs(x+50-a[0]))
        i-=i[z]
        peakys = np.append(peakys,y-i[z])
        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(401,418),_plot=False))
        plot(a[0],i)
   
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    ax1.plot(pHs)
    
    ax2.plot(peakxs,'o-')
    ax3.plot(peakys,'o-')
    ax4.plot(halfwidths,'o-')
    
    ax1.set_ylabel('pH')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
    
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    
    ax1.plot(pHs,peakys*halfwidths)
    ax2.plot(pHs,peakxs,'o-')
    ax3.plot(pHs,peakys,'o-')
    ax4.plot(pHs,halfwidths,'o-')
    
    ax1.set_ylabel('first exicton areas (relative)')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
  
 
    figure()
    for i in a[29:48:2]:
        plot(a[0],i)
    return 0    
    
def March5_etching():
    """Time dependent etching of PPA capped QDs from immediately after the PPA exchange"""
    
    a = loadtxt('160305/160305_stability measurement.csv', delimiter = ',',skiprows = 2, usecols =(0,)+tuple(range(11,75,2)),unpack = True)
    times = array([0,30,60,90,120,165,195,215])
    peakxs = array([])
    peakys = array([])
    halfwidths  = array([])

    a[1:]-=transpose([a[1:,0]])
    a[1:]=a[range(1,33,4)+range(2,33,4)+range(3,33,4)+range(4,33,4)]
    figure()
  
    for i in a[1:]:
       
        (x,y) = UVVistools.findpeak(a[0],i,(401,418))
        peakxs = np.append(peakxs,x)
        z  = argmin(abs(x+50-a[0]))
        i-=i[z]
        peakys = np.append(peakys,y-i[z])
        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(401,418),_plot=False))
        plot(a[0],i)
    peakys = peakys.reshape((4,-1))
    peakxs = peakxs.reshape((4,-1))
    halfwidths = halfwidths.reshape((4,-1))
    

   
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
#    
    for i in range(4) :
        ax1.plot(times,peakys[i]*halfwidths[i])
        ax2.plot(times,peakxs[i],'o-')
        ax3.plot(times,peakys[i],'o-')
        ax4.plot(times,halfwidths[i],'o-')
#    
    ax1.set_ylabel('first exicton areas (relative)')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')

    return 0  
    
def March7_etching():
   
    """Time dependent etching of PPA capped QDs from immediately after the PPA exchange. March 7"""
    
    a = loadtxt('160307/160307_CdSPPA-Stability.csv', delimiter = ',',skiprows = 2, usecols =(0,)+tuple(range(3,79,2)),unpack = True)
    times = array([0,30,60,90,120,180,240,300,345,405])
    peakxs = array([])
    peakys = array([])
    halfwidths  = array([])

    a[1:]-=transpose([a[1:,0]])
    b=a[array([0,
            1,5,7,11,15,19,23,27,31,35,  ### dark N2
            1,4,8,12,16,20,24,28,32,36, ### dark air
            1,6,9,13,17,21,25,29,33,37, ### light air
           # 2,3,10,14,18,22,26,30,34, ### Cd-added
            1,5,7,11,15,19,26,30,34,38])] ### light N2]
    
    figure()
  
    for i in b[1:]:
       
        (x,y) = UVVistools.findpeak(a[0],i,(390,418))
        peakxs = np.append(peakxs,x)
        z  = argmin(abs(x+50-a[0]))
        i-=i[z]
        peakys = np.append(peakys,y-i[z])
        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(390,418),_plot=False))
        plot(a[0],i)
    peakys = peakys.reshape((4,-1))
    peakxs = peakxs.reshape((4,-1))
    halfwidths = halfwidths.reshape((4,-1))
    

   
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
#    
    for i in range(3) :
        ax1.plot(times,peakys[i]*halfwidths[i])
        ax2.plot(times,peakxs[i],'o-')
        ax3.plot(times,peakys[i],'o-')
        ax4.plot(times,halfwidths[i],'o-')
    ax1.plot([0,60,120,165,225],peakys[3,5:]*halfwidths[3,5:])
    ax2.plot([0,60,120,165,225],peakxs[3,5:],'o-')
    ax3.plot([0,60,120,165,225],peakys[3,5:],'o-')
    ax4.plot([0,60,120,165,225],halfwidths[3,5:],'o-')
#    
    ax1.set_ylabel('first exicton areas (relative)')
    ax1.legend(['dark N2', 'dark air', 'light air', 'light N2'])
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
    
    figure()
    for i in b[31:]:
        d =  UVVistools.findpeak(a[0],i,(401,418))[0]
        plot(a[0]-d,i)

    return 0    
    
def March9titration():
    """pH Titration  of  CdS dots PPA capped. Samples kept in the dark during the experiment"""
    b = loadtxt('160309/160309_TitrationAdditionsAndpHs.csv',delimiter=',',skiprows = 1, unpack = True)
    pHs = b[0]
    volumes = b[1]
    volumes[8:]-=0.002
    volfix = volumes/np.roll(volumes,1)
    volfix[0]=1
    volfix[8]=1.003
    print volfix
    
   
    
    a = loadtxt('160309/160309titration.csv', delimiter = ',',skiprows = 2, usecols =(2,)+tuple(range(3,3+2*len(pHs),2)),unpack = True)
    
    
    peakxs = array([])
    peakys = array([])
    secondexciton=array([])
    thirdexciton=array([])
    halfwidths  = array([])
    a[1:]*=transpose([volfix])
    a[1:]-=transpose([a[1:,0]])
    
    
    figure()
    for i in a[1:]:
        
        (x,y) = UVVistools.findpeak(a[0],i,(410,425))
        peakxs = np.append(peakxs,x)
        z  = argmin(abs(x+50-a[0]))
        i-=i[z]
        peakys = np.append(peakys,y-i[z])
        secondexciton = np.append(secondexciton,i[z+38])
        thirdexciton = np.append(thirdexciton,i[z+68])
        halfwidths = np.append(halfwidths, UVVistools.HWHM((a[0],i),(410,425),_plot=False))
        plot(a[0],i)
   
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    ax1.plot(pHs)
    
    ax2.plot(peakxs,'o-')
    ax3.plot(peakys,'o-')
    ax3.plot(secondexciton,'o-')
    ax3.plot(thirdexciton,'o-')
    ax4.plot(halfwidths,'o-')
    
    ax1.set_ylabel('pH')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
    
    fig = figure()
    ax1 =fig.add_subplot(411)
    ax2 =fig.add_subplot(412)
    ax3 =fig.add_subplot(413)
    ax4 =fig.add_subplot(414)
    
    
    ax1.plot(pHs,peakys*halfwidths)
    ax2.plot(pHs,peakxs,'o-')
    ax3.plot(pHs,peakys,'o-')
    ax3.plot(pHs,secondexciton,'o-')
    ax3.plot(pHs,thirdexciton,'o-')
    ax4.plot(pHs,halfwidths,'o-')
    
    figure()
 
    plot(pHs,peakxs,'o-')
    plot(pHs[-10:], peakxs[-10:], 'ro-')
    
    
    
    ax1.set_ylabel('first exicton areas (relative)')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('peak absorbance')
    ax4.set_ylabel('HWHM')
    
 
    figure()
    for i in a[4:16:1]:
        d =  UVVistools.findpeak(a[0],i,(410,425))[0]
        plot(a[0]-d,i)
  
        
    figure()
    KOHmol=b[8]
    plot(halfwidths, peakxs)
    

    
    return 0
def March9_fit():
    """pH Titration  of  CdS dots PPA capped. Samples kept in the dark during the experiment"""
    b = loadtxt('160309/160309_TitrationAdditionsAndpHs.csv',delimiter=',',skiprows = 1, unpack = True)
    pHs = b[0]
    volumes = b[1]
    volumes[8:]-=0.002
    volfix = volumes/np.roll(volumes,1)
    volfix[0]=1
    volfix[8]=1.003
    print volfix
    
   
    
    a = loadtxt('160309/160309titration.csv', delimiter = ',',skiprows = 2, usecols =(2,)+tuple(range(3,3+2*len(pHs),2)),unpack = True)
    
    
    peakxs = array([])
    peakys = array([])
    peakx2=array([])
    peakx3 = array([])
    secondexciton=array([])
    thirdexciton=array([])
    halfwidths  = array([])
    
    hw1= array([])
    hw2= array([])
    hw3= array([])
    a[1:]*=transpose([volfix])
    a[1:]-=transpose([a[1:,0]])


#     
    figure()
    for i in a[1:]:
        
        (x,y) = UVVistools.findpeak(a[0],i,(410,425))
       
        wguess = 200
        Aguess = 0.4
        guess =  [Aguess,x,wguess,Aguess,x-38,wguess,Aguess,x-68,wguess,1,304,wguess,-0.001,0]#-1e-3,0]#1e-10,1e-8,
        guess =  [Aguess,Aguess,1,10,x,x-38,x-68,320,wguess,wguess,wguess,5000,]
        

        
        
        xfit = a[0][350:500]
        yfit  = i[350:500]
        z  = argmin(abs(x+50-a[0]))
        i-=i[z]
        fitres = fitspectrum(RamanSpectrum(pandas.Series(yfit,xfit)),(333,450),'xGaussianNoBase', guess)

       
       
        
        peakys = np.append(peakys,fitres.params[0][3])#-i[z])
        
        peakxs = np.append(peakxs,fitres.params[0][7])
        peakx2 = np.append(peakx2,fitres.params[0][6])
        peakx3 = np.append(peakx3,fitres.params[0][5])
        halfwidths = np.append(halfwidths,fitres.areas[3])# UVVistools.HWHM((a[0],i),(410,425),_plot=False))
        
        hw1 = np.append(hw1,sqrt(fitres.params[0][11]))
        hw2 = np.append(hw2,sqrt(fitres.params[0][10]))
        hw3 = np.append(hw3,sqrt(fitres.params[0][9]))
        secondexciton = np.append(secondexciton,fitres.areas[2])
        thirdexciton = np.append(thirdexciton,fitres.areas[1])
        
        plot(a[0],i,'--y')
        plot(fitres.x,fitres.peaks[3]+fitres.peaks[0],'k')
        plot(fitres.x,fitres.peaks[2]+fitres.peaks[0],'r')
        plot(fitres.x,fitres.peaks[1]+fitres.peaks[0],'b')

    for x in [3,27,39]:
        halfwidths[x:]/=halfwidths[x]
        secondexciton[x:]/=secondexciton[x]
        thirdexciton[x:]/=thirdexciton[x]
        
    fig = figure()
    ax1 =fig.add_subplot(141)
    ax2 =fig.add_subplot(142)
    ax3 =fig.add_subplot(143)
    ax4 =fig.add_subplot(144)
    
    ax1.plot(pHs)
    
    ax2.plot(peakxs,'o-')
    ax2.plot(peakx2,'o-')
    ax2.plot(peakx3,'o-')
    
    ax3.plot(halfwidths,'o-')
    ax3.plot(secondexciton,'o-')
    ax3.plot(thirdexciton,'o-')
    
    ax4.plot(hw1,'o-')
    ax4.plot(hw2,'o-')
    ax4.plot(hw2,'o-')
    
    
    ax1.set_ylabel('pH')
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('noramlized peak areas')
    ax4.set_ylabel('HWHM')
   
    fig = figure()
    ax2 =fig.add_subplot(131)
    ax3 =fig.add_subplot(132)
    ax4 =fig.add_subplot(133)
    
    

    ax2.plot(pHs,peakxs,'o-')
    ax2.plot(pHs,peakx2,'o-')
    ax2.plot(pHs,peakx3,'o-')
    
    ax3.plot(pHs,halfwidths/halfwidths[0],'o-')
    ax3.plot(pHs,secondexciton/secondexciton[0],'o-')
    ax3.plot(pHs,thirdexciton/thirdexciton[0],'o-')
    
    ax4.plot(pHs,hw1,'o-')
    ax4.plot(pHs,hw2,'o-')
    ax4.plot(pHs,hw3,'o-')
    
 
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('normalized peak area')
    ax4.set_ylabel('peak hwhm')
    ax4.legend(['first', 'second', 'third'])
    
    figure()
 
    plot(pHs,peakxs,'o-')

    return 0   
    
def March10titration():
    """pH Titration  of  CdS dots PPA capped. Samples kept in the dark during the experiment"""
    b = loadtxt('160310/160310_TitrationAdditionsAndpHs.csv',delimiter=',',skiprows = 1, unpack = True)
    pHs = b[0]
    volumes = b[1]
    volumes/=volumes[0]
    volfix = volumes/np.roll(volumes,1)
    
    
    
   
    
    a = loadtxt('160310/160310_PPATitration.csv', delimiter = ',',skiprows = 2, usecols =(2,)+tuple(range(3,3+2*len(pHs),2)),unpack = True)
    
    
    peakxs = array([])
    peakys = array([])
    peakx2=array([])
    peakx3 = array([])
    secondexciton=array([])
    thirdexciton=array([])
    halfwidths  = array([])
    baselines = array([])
    
    hw1= array([])
    hw2= array([])
    hw3= array([])
    
    a[1:]*=transpose([volfix])
    a[1:]-=transpose([a[1:,0]])
    
   

    figure()
    for i in a[1:]:
        
        (x,y) = UVVistools.findpeak(a[0],i,(410,425))
       
        wguess = 200
        Aguess = 0.4
        guess =  [Aguess,x,wguess,Aguess,x-38,wguess,Aguess,x-68,wguess,1,304,wguess,-0.001,0]#-1e-3,0]#1e-10,1e-8,
        guess =  [Aguess,Aguess,1,10,x,x-38,x-68,320,wguess,wguess,wguess,5000,]
        

        
        
        xfit = a[0][350:500]
        yfit  = i[350:500]
        z  = argmin(abs(x+50-a[0]))
        
        i-=i[z]
        fitres = fitspectrum(RamanSpectrum(pandas.Series(yfit,xfit)),(333,450),'xGaussianNoBase', guess)

       
       
        
        peakys = np.append(peakys,fitres.params[0][3])#-i[z])
        
        peakxs = np.append(peakxs,fitres.params[0][7])
        peakx2 = np.append(peakx2,fitres.params[0][6])
        peakx3 = np.append(peakx3,fitres.params[0][5])
        halfwidths = np.append(halfwidths,fitres.areas[3])# UVVistools.HWHM((a[0],i),(410,425),_plot=False))
        
        hw1 = np.append(hw1,sqrt(fitres.params[0][11]))
        hw2 = np.append(hw2,sqrt(fitres.params[0][10]))
        hw3 = np.append(hw3,sqrt(fitres.params[0][9]))
        secondexciton = np.append(secondexciton,fitres.areas[2])
        thirdexciton = np.append(thirdexciton,fitres.areas[1])
        x1 = argmin(abs(x+50-a[0]))
        x2 = x1-50
  
        slope = -numpy.polynomial.polynomial.polyfit(a[0,x2:x1],i[x2:x1],1)[1]
        print slope
        baselines = append(baselines, slope)
        plot(a[0],i,'--y')
        plot(fitres.x,fitres.peaks[3]+fitres.peaks[0],'k')
        plot(fitres.x,fitres.peaks[2]+fitres.peaks[0],'r')
        plot(fitres.x,fitres.peaks[1]+fitres.peaks[0],'b')

#
#    fig = figure()
#    ax1 =fig.add_subplot(411)
#    ax2 =fig.add_subplot(412)
#    ax3 =fig.add_subplot(413)
#    ax4 =fig.add_subplot(414)
#    
#    ax1.plot(pHs)
#    
#    ax2.plot(peakxs,'o-')
#    ax2.plot(peakx2,'o-')
#    ax2.plot(peakx3,'o-')
#    
#    ax3.plot(halfwidths,'o-')
#    ax3.plot(secondexciton,'o-')
#    ax3.plot(thirdexciton,'o-')
#    
#    ax4.plot(hw1,'o-')
#    ax4.plot(hw2,'o-')
#    ax4.plot(hw2,'o-')
#    
#    
#    ax1.set_ylabel('pH')
#    ax2.set_ylabel('peak center (nm)')
#    ax3.set_ylabel('noramlized peak areas')
#    ax4.set_ylabel('HWHM')
    figure()
    plot(pHs, baselines)
    fig = figure()
    ax2 =fig.add_subplot(311)
    ax3 =fig.add_subplot(312)
    ax4 =fig.add_subplot(313)
    
    

    ax2.plot(pHs,peakxs,'o-')
    ax2.plot(pHs,peakx2,'o-')
    ax2.plot(pHs,peakx3,'o-')
    
    ax3.plot(pHs,halfwidths/halfwidths[0],'o-')
    ax3.plot(pHs,secondexciton/secondexciton[0],'o-')
    ax3.plot(pHs,thirdexciton/thirdexciton[0],'o-')
    
    ax4.plot(pHs,hw1,'o-')
    ax4.plot(pHs,hw2,'o-')
    ax4.plot(pHs,hw3,'o-')
    
 
    ax2.set_ylabel('peak center (nm)')
    ax3.set_ylabel('normalized peak area')
    ax4.set_ylabel('peak hwhm')
    ax4.legend(['first', 'second', 'third'])
    
    figure()
 
    plot(pHs,peakxs,'o-')

    return 0

def CdOH2_solubility():
    K = 7.2E-15
    pH = linspace(5,14,100)
    print pH
    Cd = K/10**(-2*(14.0-pH))
    print Cd
    plot(pH,Cd ) 
    return 0