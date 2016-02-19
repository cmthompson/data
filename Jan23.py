# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 18:23:47 2015

@author: chris
"""


from ramanTools.RamanSpectrum import *


    
def Jan23():#### Work up some XPS data on CdSe deposited on Ag with varying concentration of CdSe.  
    import pandas
    pandas.options.display.mpl_style = False
    rough = loadtxt('/home/chris/Documents/DataWeiss/150123/RoughSilverData.csv',
            delimiter = ',',
            skiprows = 1,
            unpack = True)
    smooth = loadtxt('/home/chris/Documents/DataWeiss/150123/SmoothSilverData.csv',
            delimiter = ',',
            skiprows = 1,
            unpack = True)
            
    longtime = loadtxt('/home/chris/Documents/DataWeiss/150127/silver.csv',
            delimiter = ',',
            skiprows = 1,
            unpack = True)
            
    print 'first the rough sample, then smooth silver, then long time'
    print '-------------------------------------'
    for i in range(3):
        if i == 0:
            spec = rough
            title = 'rough silver'
        elif i == 1:
            spec = smooth
            title = 'smooth silver'
        elif i ==2:
            spec = longtime
            title = 'long time'
            
        
        
        figure().suptitle(title,size = 24)
        subplot(331)
        a = RamanSpectrum(pandas.Series(spec[11],spec[10]))
        a.plot()
        cadmium_area = calc_area(a,(408,401))+calc_area(a,(416,409))
        cadmium5_2_area = calc_area(a,(408,401))
        
        start = argmin(abs(array(a.index)-408))
        end = argmin(abs(array(a.index)-401))
        print start,end
        xs =a.index[start:end]
        ys= a.values[start:end]
        
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = slope*(xs-xs[0])+ys[0]
        fill_between(xs,baseline,ys)
        annotate('Cd 3d', (0.1,0.1),xycoords = 'axes fraction')
        
        xlim((390,420))
        
        
        
        subplot(332)
        a = RamanSpectrum(pandas.Series(spec[9],spec[8]))
        silver_area = calc_area(a,(378,365))#+calc_area(a,(376,372))
        
        a.plot()
        start = argmin(abs(array(a.index)-378))
        end = argmin(abs(array(a.index)-365))
        print start,end
        xs =a.index[start:end]
        ys= a.values[start:end]
        
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = slope*(xs-xs[0])+ys[0]
        fill_between(xs,baseline,ys)
        axis('on')
        xlim((350,390))
        annotate('Ag 3d', (0.1,0.1),xycoords = 'axes fraction')
        
        subplot(333)
        a = RamanSpectrum(pandas.Series(spec[3],spec[2]))
        
        a.plot()
       
        selenium_area = calc_area(a,(56,52))
        start = argmin(abs(array(a.index)-56))
        end = argmin(abs(array(a.index)-52))
        print start,end
        xs =a.index[start:end]
        ys= a.values[start:end]
        
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = slope*(xs-xs[0])+ys[0]
        fill_between(xs,baseline,ys)
        xlim((45,65))
        annotate('Se 3d', (0.1,0.1),xycoords = 'axes fraction')
        #print 'selenium atoms per silver atom:', (sulfur_area/0.6666)/(silver_area/5.987)
        #### sensitivities for 54.6 degrees,  come from Handbook of XPS, p. 253
        
        subplot(334)
        a = RamanSpectrum(pandas.Series(spec[7],spec[6]))
        
        a.plot()
        xlim((155,175))
        sulfur_area = calc_area(a,(175,155))
        annotate('S 2p', (0.1,0.1),xycoords = 'axes fraction')
        
        
        subplot(335)
        a = RamanSpectrum(pandas.Series(spec[5],spec[4]))
        
        a.plot()
        start = 0#argmin(abs(array(a.index)-135))
        end = argmin(abs(array(a.index)-130))
       
        xs =a.index[start:end]
        ys= a.values[start:end]
        
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = array([ys[-1]]*len(xs))
        phosphorous_area = sum((ys-baseline)[1:]*diff(xs))
        fill_between(xs,baseline,ys)
        xlim((120,140))
      
        
        subplot(339)
        a = RamanSpectrum(pandas.Series(spec[1],spec[0]))
        
        a.plot()
        
        
        #### sensitivities for 54.6 degrees,  come from Handbook of XPS, p. 253
    
        print 'Se atoms per silver atom', (selenium_area/0.853)/(silver_area/5.987)
        print 'Cd atoms per silver atom', (cadmium5_2_area/3.98)/(silver_area/5.987)
        print 'sulfur atoms per silver atom:', (sulfur_area/0.6666)/(silver_area/5.987)
        print 'P atoms per silver atom', (phosphorous_area/0.486)/(silver_area/5.987)
        print '__________________________________-'
    return 0
    
    
def Cd():
     
     l = linspace(0,10,50) # nm
     N_overlayer = 10 #nm^-3
     N_CdSe =18
     MFP = 1.0 #nm
     
     for i in range(0,100,20):
         N_overlayer = i
         I_Cd = (MFP*(1-exp(-l/MFP))*N_overlayer + exp(-l/MFP)*MFP*N_CdSe)
     
         I_Se = (exp(-l/MFP)*MFP*N_CdSe)
         plot(l, I_Se/I_Cd,label = str(i))
     legend()
     xlabel('Cd shell thickness')
     ylabel('Signal Se/Cd')
     return 0

def Ag():
     
     l = 4 #nm diameter of quantum dot
     
     N_CdSe =18 #Atoms Cd per nm^3
     
     MFP = 1.0 #nm
     N_Ag = 4/(0.407**3) #Atoms Ag per nm^3
     
     coverage = linspace(0,1,50)
         
     I_Cd = MFP*(1-exp(-l/MFP))*N_CdSe*coverage
 
     I_Ag = exp(-l/MFP)*MFP*N_Ag*coverage+ MFP*N_Ag*(1-coverage)
     
     plot(coverage, I_Cd/I_Ag)
     ylabel('apparent atoms Cd/atoms Ag')
     xlabel('surface coverage')
     return 0
    
def shirley(x,a,F,E):
    return cos(3.14159*a/2+(1-a)*arctan((x-E)/F))/(F**2 + (x-E)**2)**((1-a)/2)
     
   
    
    
    