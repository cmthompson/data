# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 17:37:42 2016

@author: chris
"""
import numpy
from numpy import *
from matplotlib.pyplot import *


def findpeak(x,y,rnge,_plot=False,precision=0.01):
    """Locate peak of quantum dot spectrum.  Input x, the wavelengths; and y, the absorbances; and rnge, a length-two tuple of lower and upperbound for searching for a peak.  Return length-two tuple (lambdamax and absorbance)"""
    x1 = numpy.argmin(abs(x-rnge[0]))   #### find the index of wavelength for lower bound wavelength
    x2 = numpy.argmin(abs(x-rnge[1]))   #### find the index of wavelength for upper bound wavelength
    if x1>x2:  ### make sure the indices are in the right order
        xtemp = x1
        x1=x2
        x2=xtemp
    yfit = numpy.polyfit(x[x1:x2],y[x1:x2],5)  ## Fit the y data to a 5th order polynomial
    
    xs_fit = numpy.arange(rnge[0],rnge[1],precision)  
    ys_fit = polyeval(yfit,xs_fit)
   
    if _plot:
        plot(xs_fit,ys_fit)
   
    xmax = xs_fit[argmax(ys_fit)]   ### calculate location (x-value) of the max of your polynomial fit

    ymax=max(ys_fit)   ### calculate y-value of the max of your polynomial fit
    
    return (xmax,ymax)


def QDconc(peak):
    """Calculate concentration of quantum dots from absorbance data.  Input the lambdamax and absorbance (in tuple). Give back concentration for a 1 cm cuvette"""
    if peak[0]>=500:
        print "CdSe quantum dot"
        diameter =0.0000000016122*peak[0]**4-0.0000026575*peak[0]**3+0.0016242*peak[0]**2-0.4277*peak[0]+41.57
        print 'diam', diameter
        epsilon =5857*diameter**2.65         
        print 'eps', epsilon    
        print 'conc for 1 mm cuvette', peak[1]/epsilon 
    else:  
        print "CdS quantum dot"
        diameter = -0.000000066521*peak[0]**3+0.00019557*peak[0]**2-0.092352*peak[0]+13.29
        print 'diam', diameter
        epsilon = 21536*diameter**2.3
        print 'eps', epsilon    
        print 'conc for 1 cm cuvette', peak[1]/epsilon 
        
    return peak[1]/epsilon 

def get_peak_centers():
    close('all')
    f = figure()
    ax = f.add_subplot(111)  ### create an axis for plotting the final peak values
    os.chdir('/home/chris/Desktop/Addition of Bases')  ### change to the directory containing your data
    namelist = ['1_dew_74_300eqPPA_KOH_HCl_all.csv',
                '1_dew_74_200eqPPA_KOH_HCl_all.csv']  ### put here a list of all the files you want to check out
                
                
    for name in namelist:  ### iterate over each file name

        peaks = numpy.array([])  ## initalize an array to hold the peak centers
        
       
        figure()
        print name
        a = numpy.loadtxt(name, unpack = True, skiprows = 2,delimiter=',')  ### load the csv file
       
        for i in range(0,a.shape[0],2):  ## for each set of x and y data columns....
           
            plot(a[i],a[i+1])  ###### plot the data....
            try:
                p = findpeak(a[i],a[i+1],(405,430))  ### use the 'findpeak' function to find a peak.......
            except:
                p=(0,0)  #### if you have a problem finding the peak, just give back a (0,0) result
                print 'error', name, i      ### and report the error...
            peaks = append(peaks, p[0])  ### then add the first value of your result to the array called 'peaks'
            
        
        ax.plot(peaks,label = name)  #### plot the peaks on the axis named 'ax'
        savetxt('peaks_'+name,peaks,delimiter = ',')  ### save the peak centers to a file name starting with 'peaks_'
        print peaks
    ax.legend()  #### make a legend for the axis named 'ax' 
    
    return 0
        
        
        