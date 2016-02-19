# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:00:46 2015

@author: chris
"""
from ramanTools.RamanSpectrum import *

    
def remove_CCD_interference(spectrum):  ########### remove interference noise from 3_light

    r=spectrum.copy()
    r = r.iloc[300:800]
    lambdaa = 10**7/(10**7/514.5-array(r.index))
    
    
    r.index = pandas.Float64Index(lambdaa)
    
    
  
    for i in range(3):
        r = smooth(r)
    fit = polyfit(array(r.index), r.values, 5)
    r/=polyeval(fit,array(r.index))
    xs = array(r.index)/1000
    def funct(x,A,a,b,c):return (A/1000.0)*cos(1000/(a*x**2+b*x+c))+1   #def funct(x,A,b,phase):return (1+(A)*cos(x*pi/b+phase))**2
    guess = [1.7,3,1,0]
    res = scipy.optimize.curve_fit(funct,xs,r.values,guess)
    
    
    newspectrum=spectrum.copy()
    lambdaa = 10**7/(10**7/514.5-array(newspectrum.index))
    
    fixfactor = funct(lambdaa/1000,*res[0])
    
    newspectrum.values[:]/=fixfactor

    return newspectrum
    

    

   

def darkcounts():  ### shows the digitized nature of dark spectra and very low signal spectra due to ADC gain setting.  
    figure()
    a = loadtxt('/home/chris/Documents/DataWeiss/150304/dark measurement 100 scans 1s each.txt',delimiter=',', unpack=True)
    b = loadtxt('/home/chris/Documents/DataWeiss/150304/18_ dark specrtrum 10s accum.txt',delimiter=',', unpack=True)
    c = loadtxt('/home/chris/Documents/DataWeiss/150304/19_ dark specrtrum 50s accum.txt',delimiter=',', unpack=True)
    d = loadtxt('/home/chris/Documents/DataWeiss/150304/20_ dark specrtrum 100s accum.txt',delimiter=',', unpack=True)    
    
    counts = array(list((mean(i.flatten()) for i in [a,b,c,d])))
    print counts
    plot([1,10,50,100],counts)
    title('total dark baseline for increasing accumulation times')
    xlabel('accumulation time (s)')
    ylabel('dark baseline (counts)')

    
    a = loadtxt('/home/chris/Documents/DataWeiss/150304/10_dark2.txt',delimiter=',', unpack=True)
    b = loadtxt('/home/chris/Documents/DataWeiss/150304/11_dark gain on 3.txt',delimiter=',', unpack=True)
    c = loadtxt('/home/chris/Documents/DataWeiss/150304/15_10s accums gain 1.txt',delimiter=',', unpack=True)
    d = loadtxt('/home/chris/Documents/DataWeiss/150304/16 10s accums gain 3.txt',delimiter=',', unpack=True)
    print list((i.shape for i in [a,b,c,d]))
    figure()
    plot(a[:,0])
    plot(b[:,00]/10)
    figure()
    av = mean(a,axis=1)
    subplot(421)
    bins = arange(-.105,.115,0.0025)
    hist(a[100]/mean(a[100])-1,bins=bins )
   
    print 'standard deviation 1 s acuumulations, gain=1:',std(a[100]/mean(a[100])-1)
    xlim(-.105,0.105)
    subplot(422) 
    hist(a[150]/mean(a[150])-1,bins=bins )
    print 'standard deviation 1 s acuumulations, gain=1:',std(a[150]/mean(a[150])-1)
    xlim(-.105,0.105)
    
    
    
    av = mean(b,axis=1)
    subplot(423)
    
    hist(b[100]/mean(b[100])-1,bins=bins )
    xlim(-.105,0.105)
    print 'standard deviation 1 s acuumulations, gain=3:',std(b[100]/mean(b[100])-1)
    subplot(424) 
    hist(b[150]/mean(b[150])-1,bins=bins )
    print 'standard deviation 1 s acuumulations, gain=3:',std(b[150]/mean(b[150])-1)
    xlim(-.105,0.105)
    
    av = mean(c,axis=1)
    subplot(425)
    
    hist(c[100]/mean(c[100])-1,bins=bins )
    xlim(-.105,0.105)
    print 'standard deviation 10 s acuumulations, gain=3:',std(c[100]/mean(c[100])-1)
    subplot(426) 
    hist(c[150]/mean(c[150])-1,bins=bins )
    print 'standard deviation 10 s acuumulations, gain=1:',std(c[150]/mean(c[150])-1)
    xlim(-.105,0.105)
    
    av = mean(b,axis=1)
    subplot(427)
    
    hist(d[100]/mean(d[100])-1,bins=bins )
    xlim(-.105,0.105)
    print 'standard deviation 10 s acuumulations, gain=3:',std(d[100]/mean(d[100])-1)
    subplot(428) 
    hist(d[150]/mean(d[150])-1,bins=bins )
    print 'standard deviation 10 s acuumulations, gain=3:', std(d[150]/mean(d[150])-1)
    xlim(-.105,0.105)
    
    return 0
    
    
    
    
    

    
def remove_dust(a,centers = (92,452,479,595,663,941),demo = False,blind = False):
    
    a = a.copy()
    
    def function(x,A1,w1,G1): return 1 - A1*exp(-(x-w1)**2/G1)
        
    xs = array(a.index)
    ys = a.values
    fit = polyfit(xs[570:630],ys[570:630],4)
    smoothed = polyeval(fit,xs[570:630])
    calibration_point = argmin(a.values[570:630]-smoothed)+570
    print 'calib point', calibration_point
    pix_adjust = calibration_point-595
        
#    dust_peak_parameters = ([  -0.07,   87 + pix_adjust,  5],
#                            [  -1.21484283e-02,  451+ pix_adjust, 4],
#                            [  -1.13809633e-02,  479+ pix_adjust,  2],
#                            [  -1.59194998e-02,  595+ pix_adjust,   6],
#                            [  -7.46384750e-03,   623+ pix_adjust,   3],
#                            [  -7.46384750e-03,   665+ pix_adjust,   3.4],
#                            [  -0.01,   940+ pix_adjust, 5])   #### relative amplitue, center (in pixel), width of each peak. 
      
      
    dust_peak_parameters =( [ -1.33573613e-02,  87 + pix_adjust,  10],
            [ -1.07265758e-02,   452 + pix_adjust,   5],
            [  1.04830060e-02,  484 + pix_adjust ,  2],
            [ -1.69072443e-02,   595 + pix_adjust,  6],
            [ -1.47284158e-02,   623 + pix_adjust,  3],
            [ -1.71921489e-02,   665 + pix_adjust,   3],
            [ -1.82566862e-02 ,  940+ pix_adjust,   10])
    if blind == True:
        

        for z in dust_peak_parameters:
            w =round(z[2])*3
            center = z[1]
            x = arange(center-w,center+w,1)
            a.iloc[center-w:center+w]*=function(x,*z)  
        
        
    else:        
        for z in dust_peak_parameters:
            w =round(z[2])*3
            center = z[1]
            x = arange(center-w,center+w,1)
            y = a.values[center-w:center+w]
            
            slope =(y[-1]-y[0])/(x[-1]-x[0])
            b = y[0]-slope*x[0]
            baseline = slope*(x)+b
            
            
         
            z[0] = y[-1]*z[0]
           
            try:
                result = scipy.optimize.curve_fit(function,x,baseline/y,z)
                
                
                
                print result[0]
            except RuntimeError:
                print('Fit Failed')
                continue
            z = list(result[0])
            
            a.iloc[center-w:center+w]*=function(x,*z)
    if demo:
            print z
            plot(x,function(x,*z))
    return a

def removecorrelatednoise(spectrum):
    
    spectrum = spectrum._copy()
    b=loadtxt('/home/chris/Documents/DataWeiss/CCD Pixel-to-Pixel Correction Factor.txt')
    if b.size != spectrum.values.size:
        raise ValueError, "sizes are wrong"
        return spectrum
    b-=1    
    xs = array(spectrum.index)
    ys= spectrum.values
  
    average = SGsmooth(xs,ys)
    
    noise=(average)/(ys)-1
    plot(noise*1000)
    plot(b*1000)
    
    

    noise=(average)/(ys)-1
    corr_factor =  correlate(b,noise)/sqrt(sum(noise**2)*sum(b**2))
    print 'before:', corr_factor
        
    spectrum*=(1+b)
   
    
    
    plot(noise*1000)
    
    corr_factor =  correlate(b,noise)/sqrt(sum(noise**2)*sum(b**2))
    print 'after:', corr_factor
    return spectrum
        
    
def remove_dark_baseline(spectrum):
    spectrum = spectrum._copy()
    sub = spectrum.num_frames*min(1320,spectrum.accum_time*118)
    print 'subtracted', sub 
    spectrum-=sub
    return spectrum
    
def tr():  ### test total removal procedure
    clf()
    ax1= subplot(221)
    ax2= subplot(222)
    ax3= subplot(223)
    ax4=  subplot(224)
    ax1.set_title('subtract dark')
    ax2.set_title('remove pixel noise')
    ax3.set_title('remove dust')
    ax4.set_title('remove waves')
    
    names  = ['/home/chris/Documents/DataWeiss/150304/3_ light.SPE',
              '/home/chris/Documents/DataWeiss/150304/4_light.SPE',
              '/home/chris/Documents/DataWeiss/150228/RhB 500sec full power_filter.SPE',
              '/home/chris/Documents/DataWeiss/150227/10 Rhb 1100 grating.SPE',
              '/home/chris/Documents/DataWeiss/150304/7_1800 grating.SPE',
              '/home/chris/Documents/DataWeiss/150304/8.SPE']
    names = ['/home/chris/Documents/DataWeiss/150304/9.SPE' ]
    for name in names:
        e = RamanSpectrum(name)
        if e is None:
            continue
        print e.num_frames, e.accum_time
        
        
        a = remove_dark_baseline(e)
        
        
        ax1.plot(e.values,color = 'b')
        ax1.plot(a.values,color = 'r')
            
        
        b = removecorrelatednoise(a)    
        
        
        ax2.plot(a.values,color = 'b')
        ax2.plot(b.values,color = 'r')
        
        c= remove_dust(b,blind=True,demo=False)
    
      
        ax3.plot(b.values,color = 'b')
        ax3.plot(c.values,color = 'r')
        
        d= remove_CCD_interference(c)
        
       
        ax4.plot(c.values,color = 'b')
        ax4.plot(d.values,color = 'r')
      
    
    
    return 0
    
def c():
    width =50
    name = '/home/chris/Documents/DataWeiss/150304/9.SPE' 
    input_array = RamanSpectrum(name)
    input_array = removecorrelatednoise(input_array)
    plot(input_array)
    return input_array
    
    