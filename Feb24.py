# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 17:11:22 2015

@author: chris
"""
import imp
import pandas
simpleFresnel = imp.load_source('respace_x', '/home/chris/PyScripts/Old stuff/simpleFresnel.py')
from pandas import Series
from SFGMeFiles.SFG_Analysis import SG_Smooth
import pdb
from RamanTools3 import *

NA = imp.load_source('remove_dust', '/home/chris/Documents/Subgroup Stuff/NoiseAnalysis_2.py')

def Lamp(): #calibration of the lamp
     global correctionVis1
     os.chdir('/home/chris/Documents/DataWeiss/150224')
     Vis1 = loadtxt('VIs1', unpack=True,skiprows = 0,delimiter='\t')
     Vis2 = loadtxt('Vi2', unpack=True,skiprows = 0,delimiter='\t')
     Vis3 = loadtxt('VIs3', unpack=True,skiprows = 0,delimiter=',')
     IR1 = loadtxt('NIR1', unpack=True,skiprows = 0,delimiter='\t')
     IR2 = loadtxt('NIR2', unpack=True,skiprows = 0,delimiter='\t')
     lamp = loadtxt('Lamp.csv', unpack=True,skiprows = 1,delimiter=',')
     
     Vis1[1]/=max(Vis1[1])
     Vis2[1]/=max(Vis2[1])
     Vis3[1]/=max(Vis3[1])
     IR1[1]/=max(IR1[1])
     IR2[1]/=max(IR2[1])
     
    #lamp_s = Series(lamp[1],lamp[0])
    # lamp_s = lamp_s.reindex(arange(200,1500,1),method = 'ffill')
     
     (lampnew_x,lampnew_y) = simpleFresnel.respace_x(lamp[0],lamp[1],arange(280,1500,1),_plot = False)
     #Vis1[1] = SGsmooth(Vis1[0],Vis1[1])
     window_len=5
     for i in range(15):
         s=r_[Vis1[1,window_len-1:0:-1],Vis1[1],Vis1[1,-1:-window_len:-1]]
         w=ones(window_len,'d')
         Vis1[1,:] =convolve(w/w.sum(),s,mode='valid')[(window_len-1)/2:-(window_len-1)/2]
         
         
     for i in range(15):
         s=r_[IR1[1,window_len-1:0:-1],IR1[1],IR1[1,-1:-window_len:-1]]
         w=ones(window_len,'d')
         IR1[1,:] =convolve(w/w.sum(),s,mode='valid')[(window_len-1)/2:-(window_len-1)/2]
         
     plot(Vis1[0],Vis1[1],label = 'vis1') 
     plot(IR1[0],IR1[1],label='IR1')
     #plot(IR2[0],IR2[1])
     plot(lampnew_x,lampnew_y,label='lamp')
    ##################
     x = argmin(abs(280-lampnew_x))
     y = argmin(abs(1000-lampnew_x))
     print x,y
     correctionVis = Vis1[1,0:721]/lampnew_y[x:y+1]
     correctionVis/=max(correctionVis)
     Viscorr = polyfit(arange(711),correctionVis[10:721],6)
     correctionVis1 = pandas.Series(polyeval(Viscorr, arange(711)),Vis1[0,10:721])
     
     correctionVis1.plot(marker = '.',label = 'correctionVis')
     ####################################
     x = argmin(abs(min(IR1[0])-lampnew_x))
     y = argmin(abs(max(IR1[0])-lampnew_x))
     print lampnew_x[x], IR1[0,0]
     print lampnew_x[y], IR1[0,-2]
     
     
     
     
     correctionIR = IR1[1,0:-1]/lampnew_y[x:y+1]
     correctionIR/=max(correctionIR)
     print arange(IR1[0].size-11).size
     print correctionIR[10:].size
    
     IRcorr = polyfit(arange(IR1[0].size-11),correctionIR[10:],6)
     
     correctionIR1 = pandas.Series(polyeval(IRcorr, arange(IR1[0].size-11)),IR1[0,10:-1])
     correctionIR1.plot(marker = '.',label = 'correctionIR')
     
     
     legend()
     
     
     
     s= loadtxt('/home/chris/Documents/DataWeiss/150225/Spectra of RhodamineB and CresylViolet.csv',
                  skiprows=1,
                  unpack=True,
                  delimiter=',')
     wl = s[0]
     RBbase = s[1]/max(s[1])
     RB=s[2]/max(s[2])
     CVbase = s[3]/max(s[3])
     CV = s[4]/max(s[4])
     
#     figure()
#     
#     correctionVis1.plot()
#     plot(wl, 1/correctionVis1[500:700])
     
     RBref = pandas.Series(array([17.55, 32.7,54.69,77.96, 95.02, 99.85, 92.86,79.99,65.65,52.92,42.59,35.4,30.63,28,26.08])/100,
            array([545, 550,555,560,565,570,575,580,585,590,595,600,605,610,615]))
     RBref2 =loadtxt('/home/chris/Documents/Literature/Fluorescence standards/RhodamineB',
                     delimiter = '\t',
                     comments = '#',
                     unpack = True)
     RBref2 = pandas.Series(RBref2[1]/max(RBref2[1]), RBref2[0])
     RBcorrected = RB/correctionVis1[500:700]
     RBcorrected/=max(RBcorrected)
     RBcorrected.to_csv('/home/chris/Documents/DataWeiss/RhodamineB.csv')
     
     CVcorrected = CV/correctionVis1[500:700]
     CVcorrected/=max(CVcorrected)
     
     figure()
     
     plot(wl, correctionVis1[500:700],label = 'corrfactor')
     plot(wl, RB, 'k',label='RB')
     plot(wl, RBcorrected, 'r',label = 'RBcorr')
     
     RBref.plot(marker='s', label = 'RBLit')
     
     RBref2.plot(label='RBLit2')
     legend()
     
     figure()
     plot(wl, correctionVis1[500:700],label = 'corrfactor')
     CVref =loadtxt('/home/chris/Documents/Literature/Fluorescence standards/CresylViolet',
                     delimiter = '\t',
                     comments = '#',
                     unpack = True)
     CVref = pandas.Series(CVref[1]/max(CVref[1]), CVref[0])
     plot(wl, CV, 'k', label = 'CV')
     plot(wl, CVcorrected, 'r-', label = 'CVcorr')
     CVref.plot(label = 'CVlit')
     legend()

             
     
     return 0

def calibration():
    global correctionVis1
    os.chdir('/home/chris/Documents/DataWeiss/150227')
    spec1=RamanSpectrum('5 Rhb.SPE')
    names = ['1 Rhb.SPE',
             #'2 Rhb.SPE',
             #'3 Rhb.SPE',
             #'4 Rhb.SPE',
              '5 Rhb.SPE',
               '6 Rhb.SPE']
#            #    '7 RhB.SPE'
#    #'8 Rhb 1800.SPE',
#    # '9RhB 1800.SPE'
#    '10 Rhb 1100 grating.SPE',
#             '11_Rhb 1100.SPE',
#              '12_.SPE',
#               '13_rhb high conc.SPE',
#               '14.SPE',
#               '15.SPE']
                
    clf()
               
    sum_array = zeros((1024,1))         
    ax1=subplot(221)
    ax2=subplot(222)
    ax3=subplot(223)
    r = spec1.size-1
    
    
    darksignal =0# mean(RamanSpectrum( 'dark 50 s.SPE',))*10
    spec1-=darksignal
 
    xs = array(spec1.index)
    ys= spec1.values
  
    average = SGsmooth(xs,ys)
    fit = polyfit(xs,ys,6)
    dust = polyeval(fit,xs)
    noise1=transpose([ys/average])
    dustnoise = dust/average
    fullnoise = dustnoise*noise1.flatten()
    
    ax1.plot(fullnoise)
   
    sumnoise = sum(noise1**2)
    print 'dark signal cps', darksignal/500
    
    print darksignal
    for name in names:
        
        spectrum = RamanSpectrum(name)
        def match(x,A,b):
            
            return A*x-b
        x0=[10,1000]
        res = scipy.optimize.curve_fit(match,spectrum.values,spec1.values,x0)
        darksignal = res[0][1]
        mult = res[0][0]
        
        
                
        xs = array(spectrum.index)
        ys= match(spectrum.values,mult,darksignal)
        average = SGsmooth(xs, ys)
        fit = polyfit(xs,ys,6)
        dust = polyeval(fit,xs)
       
        dustnoise = dust/average
        noise = transpose([ys/(average)])
        fullnoise = (noise.flatten())*(dustnoise)
        fullnoise = SGsmooth(xs, fullnoise)
       
        
        print name,darksignal,mult, correlate(noise[:,0]-1, noise1[:,0]-1)/sqrt(sum((noise-1)**2)*sumnoise)
        
        sum_array = append(sum_array,noise,axis=1)
        
        ax1.plot(fullnoise)
        ax2.plot(noise)
        ax3.plot(xs, ys)
        #ax3.plot(xs,average)
   
    ax3.legend(['a','b','c','d','e','f'])
    ax2.legend(list(x[-14:-9] for x in names))
    subplot(224)
    
    sum_array = sum_array[:,1:]
    CCDcorrectionfactor = 1+mean(sum_array,axis=1)
    
    errorbar(range(1024), CCDcorrectionfactor, yerr=std(sum_array,axis=1)/sqrt(len(names)))

    
    
def calibration2():
    global correctionVis1
    os.chdir('/home/chris/Documents/DataWeiss/150210')
    spec1=RamanSpectrum('1_under N2 0 min.SPE')
    names = ['3_under N2 2 min.SPE',
               '4_under N2 3 min.SPE',
               '5_under N2 4 min.SPE',
               '6_under N2 5 min.SPE',
               '7_under N2 6 min.SPE',
               '8_under N2 7 min.SPE',
               '9_under N2 8min.SPE']
    
               
    clf()
    spec1 = remove_dust(spec1)
    sum_array = zeros((1024,1))         
    ax1=subplot(221)
    ax2=subplot(222)
    ax3=subplot(223)
    r = spec1.size-1
    
    
    darksignal =0# mean(RamanSpectrum( 'dark 50 s.SPE',))*10
    spec1-=darksignal
 
    xs = array(spec1.index)
    ys= spec1.values
  
    average = SGsmooth(xs,ys)
    fit = polyfit(xs,ys,6)
    dust = polyeval(fit,xs)
    noise1=transpose([ys/average])
    dustnoise = dust/average
    fullnoise = dustnoise*noise1.flatten()
    
    ax1.plot(fullnoise)
    
    sumnoise = sum(noise1**2)
    print 'dark signal cps', darksignal/500
    
    print darksignal
    for name in names:
        
        spectrum = RamanSpectrum(name)
        def match(x,A,b):
            
            return A*x-b
        x0=[10,1000]
        res = scipy.optimize.curve_fit(match,spectrum.values,spec1.values,x0)
        darksignal = res[0][1]
        mult = res[0][0]
        
        
                
        xs = array(spectrum.index)
        ys= match(spectrum.values,mult,darksignal)
        average = SGsmooth(xs, ys)
        fit = polyfit(xs,ys,6)
        dust = polyeval(fit,xs)
       
        dustnoise = dust/average
        noise = transpose([ys/(average)])
        fullnoise = (noise.flatten())*(dustnoise)
       
        
        print name,darksignal,mult, correlate(noise[:,0]-1, noise1[:,0]-1)/sqrt(sum((noise-1)**2)*sumnoise)
        
        sum_array = append(sum_array,noise,axis=1)
        ax1.plot(dustnoise)
        ax1.plot(fullnoise)
        ax2.plot(noise)
        ax3.plot(xs, ys)
        #ax3.plot(xs,average)
   
    ax3.legend(['a','b','c','d','e','f'])
    ax2.legend(list(x[-14:-9] for x in names))
    subplot(224)
    
    sum_array = sum_array[:,1:]
    CCDcorrectionfactor = 1+mean(sum_array,axis=1)
    
    errorbar(range(1024), CCDcorrectionfactor, yerr=std(sum_array,axis=1)/sqrt(len(names)))

    
    return 0
    
     
     
def SERS():
    #pdb.set_trace()
    ref = loadtxt('/home/chris/PyScripts/SilverVisRefInd.csv', delimiter = ',', usecols = (0,1,2), skiprows = 1, unpack=True)
    l= 1240/ref[0]
    j =  complex(0,1)
    n2 = ref[1]+j*ref[2]
    e = n2**2
    e0=1.77
    g = (e-e0)/(e+2*e0)
    alpha_zz = abs(1+2*g)**4
    alpha_xz = 2*abs(1+2*g)**2*abs(1-g)**2
    alpha_xx = 4*abs(1-g)**4
    plot(l,alpha_zz,'r',label='zz')
    plot(l,alpha_xz, 'b',label='xz')
    plot(l,alpha_xx,'k',label = 'xx')
    yscale('log')
    legend()
    return 0
    
def calibration3(save_it= False):
    global correctionVis1
    ax1=subplot(221)
    ax2=subplot(222)
    ax3=subplot(223)
    os.chdir('/home/chris/Documents/DataWeiss/150228')
    spec1=RamanSpectrum('RhB 500sec full power_filter.SPE')-1320
    spec1.plot(axes =ax3)
    names = ['/home/chris/Documents/DataWeiss/150227/1 Rhb.SPE',
             '/home/chris/Documents/DataWeiss/150227/10 Rhb 1100 grating.SPE',
             'dark 50 s.SPE',
             'RhB 500sec 0_01_filter.SPE',#,
             'RhB 500sec 0_1_filter.SPE',
             'RhB 500sec full power_filter.SPE',
             '/home/chris/Documents/DataWeiss/150227/1 Rhb.SPE']
    clf()
    #spec1 = NA.remove_dust(spec1,blind=True)           
    sum_array = zeros((1024,1))         
    ax1=subplot(221)
    ax2=subplot(222)
    ax3=subplot(223)
    r = spec1.size-1
    
    
    darksignal = 500*12
    spec1-=darksignal-1
 
    xs = array(spec1.index)
    ys= spec1.values
  
    average = SGsmooth(xs,ys)
    fit = polyfit(xs,ys,6)
    dust = polyeval(fit,xs)
    noise1=transpose([(average)/(ys)])
    dustnoise = dust/average
    if save_it == True:
        savetxt('/home/chris/Documents/DataWeiss/CCD Pixel-to-Pixel Correction Factor.txt', noise1)
       
    
    fullnoise = dustnoise*noise1.flatten()
   
   
    sumnoise = sum((noise1-1)**2)
    
    
    for name in names:
        
        spectrum = RamanSpectrum(name)-1320
        spectrum = spectrum.reindex(spec1.index,fill='backfill')
       
        xs = array(spectrum.index)
        ys= spectrum.values
        average = SGsmooth(xs, ys)
        noise = transpose([(average)/(ys)])
        fullnoise = (noise.flatten())*(dustnoise)
        print name,correlate(noise[:,0]-1, noise1[:,0]-1)/sqrt(sum((noise-1)**2)*sumnoise)
        
        sum_array = append(sum_array,noise,axis=1)
        
    
        ax2.plot(noise)
        ax3.plot(xs, ys)
        #ax3.plot(xs,average)
   
    ax3.legend(['a','b','c','d','e','f'])
    ax2.legend(list(x[-14:-9] for x in names))
    ax1.legend(list(x[-14:-9] for x in names))
    subplot(224)
    
    sum_array = sum_array[:,1:]
    CCDcorrectionfactor = 1+mean(sum_array,axis=1)
    
    errorbar(range(1024), CCDcorrectionfactor, yerr=std(sum_array,axis=1)/sqrt(len(names)))

    
    return 0
    

def RhodBonRaman():
    
    
    os.chdir('/home/chris/Documents/DataWeiss/150228')
    RB1=RamanSpectrum('RhB 500sec 0_01_filter.SPE')
    RB2= RamanSpectrum('RhB 500sec 0_1_filter.SPE')
    RB3=RamanSpectrum('RhB 500sec full power_filter.SPE')
    RBref = RamanSpectrum(pandas.Series.from_csv('/home/chris/Documents/DataWeiss/RhodamineB.csv'))
    RBref.index=pandas.Float64Index(10**7/514.5-10**7/array(RBref.index))
    dark = mean(RamanSpectrum('dark 50 s.SPE'))*10
    
    
    

    RB1/=max(RB1)
    RB1.plot()
    RBref.plot()
    return RBref
    