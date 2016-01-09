# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 13:51:51 2015

@author: chris
"""
from RamanTools3 import *

import os



def calculate_enchancement():
    global r,s
    #### Calc SERS enhancement on
    os.chdir('/home/chris/Documents/DataWeiss/150109')

    r = RamanSpectrum('1_MeTOP roughened Ag_1.txt')
    s = RamanSpectrum('2_MeOTP smooth silver_1.txt')
    r.autobaseline((500,1750),order=0)
    s.autobaseline((500,1750),order=0)
    on_roughened = calc_area(r,(1050,1130))*100 #### multiply, because used filter 0.01
    on_smooth = calc_area(s,(1050,1130))
    print 'hormalized area on roughened substrate =', on_roughened  
    
    print 'hormalized area on roughened substrate =', on_smooth
    r.plot(color = 'r')
    s.plot(color = 'k')
    legend(['roughened', 'smooth (x100)'],loc=2)
    annotate('x100', (0.6,0.7), xycoords = 'axes fraction', size = 24,color = 'k')
    print 'Approximate surface enhancement =',on_roughened/on_smooth
    xlabel('Raman Shift cm$^{-1}$')
    ylabel('Intensity a.u.')
    return 0
    
def figure1():
    os.chdir('/home/chris/Documents/DataWeiss/150113')

    r = RamanSpectrum('15_CdSeMTP dropcast_1.txt')
    s = RamanSpectrum('18_1.txt')
    r = smooth(r)
    s= smooth(s)
    r+=100
    s+=100
    r.plot()
    s.plot()
    
    figure()
    r = autobaseline(r,(467,1463),order=3)+200
    r.plot()
    
    s = autobaseline(s,(467,1463),order=3)+400
    s.plot()
    
    MBT = RamanSpectrum('/home/chris/Documents/DataWeiss/141014/4_methoxythiophenol_1.csv')
    MBT-=min(MBT[0:2000])
    MBT/=max(MBT[0:2000])/500
    MBT.plot()
    
    legend(['Ag/Hexanethiol', 'Ag/Hexanethiol + CdSeMTP', 'CdMTP ref'])
    
    
    return 0
    

    
def concentrationdependence():  ### determine best conc of dots to add to silver to get signal/fluorescend. 

    os.chdir('/home/chris/Documents/DataWeiss/150114')
    
    two_x = RamanSpectrum('1_concentration 2x -highest conc_1.txt')
    two_x+=1000
    one_x = RamanSpectrum('2_1x conc_1.txt')
    one_x+=500
    _25x = RamanSpectrum('3_0_25xconc_1.txt')
    _25x*=10
    _25x-=1000
    _0625x = RamanSpectrum('5_0_0625x conc_1.txt')
    _0625x*=10
    _0625x-=1500
    
    
    two_x.plot()
    one_x.plot()
    _25x.plot()
    _0625x.plot()
    
    legend(['2x','1x','0.25x*10','0.0625x * 10'])
    xlabel('Raman Shift cm$^{-1}$')
    ylabel('Intensity a.u.')
    return 0
def processxymap():
    os.chdir('/home/chris/Documents/DataWeiss/150114')
    averagespectrum = RamanSpectrum('10_maybemap_1.txt')
    
    phononarea = array([])
    fluorescenceat1600=array([])
    fullspectrum = ndarray((1000,))
    for f in os.listdir('.'):
        if '10_maybemap' in f:
            if 'SPE' in f:
                continue
            elif f == 'maybemap_1.txt':
                
                continue
            elif f == 'maybemap.txt':
                
                continue
            
            else:
                r = RamanSpectrum(f)
                phononarea = append(phononarea,r.calc_area((200,230)))
                fluorescenceat1600 = append(fluorescenceat1600,r.values[-1]-min(r.values))
                averagespectrum+=r
    averagespectrum = _smooth(averagespectrum)
    figure()
    subplot(221)
    
    
    hist(fluorescenceat1600,bins=range(0,200,20))
    hist(phononarea,bins=range(0,200,20),color='r')
    xticks(range(0,200,20))
    subplot(223)
    averagespectrum.plot()
                
                
    return phononarea
    
def XPSworkup():#### Work up some XPS data on CdSe deposited on Ag with varying concentration of CdSe.  
    import pandas
    pandas.options.display.mpl_style = None
    conc2x = loadtxt('/home/chris/Documents/DataWeiss/150116/2xsurvey.txt',
            delimiter = ',',
            unpack = True)
    concquartx = loadtxt('/home/chris/Documents/DataWeiss/150116/quarterx.csv',
            delimiter = ',',
            unpack = True)
            
            
    a = RamanSpectrum(pandas.Series(concquartx[1],concquartx[0]))
    a.plot(label='0.25x')
    print 'Cd(5/2)/Ag area 0.25x', calc_area(a,(408,402))/calc_area(a,(380,364))
 
    a = RamanSpectrum(pandas.Series(conc2x[1],conc2x[0]))
    a.plot(label='2x')
    
    print 'Cd(5/2)/Ag area 2x', calc_area(a,(408,402))/calc_area(a,(380,364))
    
    legend()
    annotate('Cd 3d', (456,183000))
    annotate('Ag 3d', (360,381000))  
    annotate('C 1s', (288,117000)) 
    
    
    
    figure()
    subplot(221)
    a = RamanSpectrum(pandas.Series(concquartx[3],concquartx[2]))
    a.plot()
    start = argmin(abs(array(a.index)-408))
    end = argmin(abs(array(a.index)-401))
    print start,end
    xs =a.index[start:end]
    ys= a.values[start:end]
    
    slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
    baseline = slope*(xs-xs[0])+ys[0]
    fill_between(xs,baseline,ys)
    cadmium5_2_area= calc_area(a,(408,402))
    
    xlim((390,420))
    
    subplot(222)
    a = RamanSpectrum(pandas.Series(concquartx[5],concquartx[4]))
    a.plot()

    xlim((275,298))
    annotate('C 1s', (0.1,0.1),xycoords = 'axes fraction')

    subplot(223)
    a = RamanSpectrum(pandas.Series(concquartx[7],concquartx[6]))
    silver_area = calc_area(a,(378,365))#calc_area(a,(370,365))+calc_area(a,(376,372))
    
    a.plot()
    start = argmin(abs(array(a.index)-378))
    end = argmin(abs(array(a.index)-365))
    xs =a.index[start:end]
    ys= a.values[start:end]
    
    slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
    baseline = slope*(xs-xs[0])+ys[0]
    fill_between(xs,baseline,ys)
    
    
    xlim((350,390))
    annotate('Ag 3d', (0.1,0.1),xycoords = 'axes fraction')
    
    
    subplot(224)
    a = RamanSpectrum(pandas.Series(concquartx[9],concquartx[8]))
    a.plot()
    xlim((155,175))
    sulfur_area = calc_area(a,(175,155))
    annotate('S 2p', (0.1,0.1),xycoords = 'axes fraction')
   
    
    
    #print 'Se atoms per silver atom', (selenium_area/0.853)/(silver_area/5.987)
    print 'Cd atoms per silver atom', (cadmium5_2_area/3.98)/(silver_area/5.987)
    print 'sulfur atoms per silver atom:', (sulfur_area/0.6666)/(silver_area/5.987)
    
    print '__________________________________-'
    #### sensitivities for 54.6 degrees,  come from Handbook of XPS, p. 253
    return 0
    
def Jan18():### data from trying out the cryostat
    os.chdir('/home/chris/Documents/DataWeiss/150118')
    green = RamanSpectrum('11_phonon long scan.SPE')
    green.plot()
    
    red = RamanSpectrum('12b.SPE')
    s = argmin(abs(red.index-1468))
    print s
    print red.iloc[s-1:s+2]
    red.iloc[s+1:] -= red.iloc[s+1]-red.iloc[s]
   # red = _smooth(red)
    red = autobaseline(red,(131,500),order = 0)
    red = autobaseline(red,(500,950),order= 4)
    red = autobaseline(red,(950,1520),order= 1)
    red.plot()
    
    
    red = RamanSpectrum('13.SPE')
    s = argmin(abs(red.index-942.3))
  
    red.iloc[s+1:] -= red.iloc[s+1]-red.iloc[s]
    #red = _smooth(red)
    red = autobaseline(red,(525,1612),order = 0)
    #red = autobaseline(red,(500,950),order= 4)
    #red = autobaseline(red,(950,1520),order= 1)
    red.plot()
    return 0
    
    
def Jan22():### data from trying out the cryostat
    from scipy import optimize
    os.chdir('/home/chris/Documents/DataWeiss/150121')
    
    
    red = RamanSpectrum('17.SPE')
    red = normalize(red,(0,1600))
    red.plot(label='647 nm')
    
    green = RamanSpectrum('14.SPE')
    green=normalize(green,(0,1600))
    green.plot(label = '514 nm')
    
    
    cdRef=normalize(CdMeOTPRef,(0,1600))
    cdRef.plot(label = 'CdMeOTP Reference')
    
    mtpref = normalize(MeOTPRef,(0,1600))
    mtpref.plot(label = 'MeOTP Reference')
    xlim(1050,1150)
    
    legend()
    
    def singlegauss(x, A1, x1, c1):return A1*exp(-(x-x1)**2/(2*c1**2)) 
    def doublegauss(x, A1, x1, c1,A2,x2,c2):return A1*exp(-(x-x1)**2/(2*c1**2)) +A2*exp(-(x-x2)**2/(2*c2**2))
    green = autobaseline(green, (1050,1150), order = 0)
    print argmin(abs(green.index-1100))
    x = array(green.index[620:680])
    y = array(green.values[620:680])  
    print x
    
    
    r = list(optimize.curve_fit(singlegauss,x,y,[1,1087,20])[0])
    figure()
    plot(x,y,'s')
    plot(x,singlegauss(x,*r),'k')
    print r
    
    r = list(optimize.curve_fit(doublegauss,x,y,[1,1080,10,1,1090,10])[0])
    
    plot(x,doublegauss(x,*r),'r')
    
    
    plot(x,singlegauss(x,r[0],r[1],r[2]),'k.')
    plot(x,singlegauss(x,*r[3:6]),'k.')
    
    
    print r
    
    
         
    xlabel('Raman Shift cm$^{-1}$')
    ylabel('Intensity a.u.')
    return 0
    
    
    
    
