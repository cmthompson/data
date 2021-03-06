# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 09:42:41 2015

@author: chris
"""
from ramanTools.RamanSpectrum import *
from numpy import *
import numpy as np

from matplotlib.pyplot import *


def findpeak(x,y,rnge,_plot=False,precision=0.01):
    """Locate peak of quantum dot spectrum.  Input spectrum and approximate range for peak.  Return lambdamax and absorbance"""
    
    x1 = np.argmin(abs(x-rnge[0]))
    x2 = np.argmin(abs(x-rnge[1]))
    if x1>x2:
        xtemp = x1
        x1=x2
        x2=xtemp
    yfit = np.polynomial.polynomial.polyfit(x[x1:x2],y[x1:x2],5)
    
    xs_fit = np.arange(rnge[0],rnge[1],precision)
    ys_fit = np.polynomial.polynomial.polyval(xs_fit,yfit)
   
    if _plot:
        plot(xs_fit,ys_fit)
   
    xmax = xs_fit[argmax(ys_fit)]

    ymax=max(ys_fit)
    
    return (xmax,ymax)

def HWHM(z,rnge,_plot=False):
    """Calculate half-width at half max of UVVis spectrum of QDs.  Input spectrum and approximate range for peak.  Return lambdamax and absorbance"""
    precision=0.1
    x=z[0]
    y=z[1]
    x1 = np.argmin(abs(x-rnge[0]))
    x2 = np.argmin(abs(x-rnge[1]))
    if x1>x2:
        xtemp = x1
        x1=x2
        x2=xtemp
    yfit = np.polynomial.polynomial.polyfit(x[x1:x2],y[x1:x2],5)
    xs_fit = np.arange(rnge[0],rnge[1],precision)
    ys_fit = np.polynomial.polynomial.polyval(xs_fit,yfit)
   
   
   
    xmax = xs_fit[argmax(ys_fit)]

    ymax=max(ys_fit)
#    halfmax = ymax/2
#    halfwidth = abs(xs_fit[argmin(abs(ys_fit-halfmax))] - xmax)  

    
    def gaussianfunction(x,A,x0,G,m,b):return m*x+b+A*exp(-(x-x0)**2/G)
    r = scipy.optimize.curve_fit(gaussianfunction, x[0:argmin(abs(x-xmax))+4], y[0:argmin(abs(x-xmax))+4], [ymax,xmax,10,0,0])[0]
    
    halfwidth = 2*sqrt(log(2)*r[2])/2
    if _plot:
        plot(x[0:argmin(abs(x-xmax))+4], gaussianfunction (x[0:argmin(abs(x-xmax))+4],*r))
    return halfwidth
    

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
    

    
def indivQY(UVVisfile, UVViscolumn, anthracenecolumn, fluorescencefile, anthracenefluorescencefile,
            subtractfluorfile= None,
            UVVisplot=None, fluorplot=None,
            fluorescencerange = (410,473),
            excitationwavelength = 350,
            nliq=1.333,
            day=0, label=None,color = 'k',
            _plot_standard = False,
            subtract_smooth_background_for_anthracene = False
            ):
    print '-------------------------------------'
    print 'calculating fluorescence yield for', label,'file', fluorescencefile
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    if len(UVViscolumn)==1:
        numuvviscolumn = alphabet.find(UVViscolumn)
    elif len(UVViscolumn)==2:
        numuvviscolumn = alphabet.find(UVViscolumn[0])*26+alphabet.find(UVViscolumn[1])
    
        
    
    a = loadtxt(UVVisfile, delimiter = ',', unpack = True, skiprows = 1,usecols=(0,numuvviscolumn,alphabet.find(anthracenecolumn)))
    
    a[1:]-=transpose([a[1:,0]])
    
   
    anthracene=RamanSpectrum(pandas.Series(a[2][::-1],a[0][::-1]))
    dot=RamanSpectrum(pandas.Series(a[1][::-1],a[0][::-1]))
    
    if subtract_smooth_background_for_anthracene :
        anthracene.smoothbaseline((290,300),(390,400))
    
    anthracene[:]-=anthracene[389]
    
    anthraceneabsorbance350= anthracene[excitationwavelength]
    absvalues = dot[excitationwavelength]
    
    nE = 1.359
    nQ = 1.44
    nW = 1.333  ## refractive index water    
 
    
    
    a = loadtxt(anthracenefluorescencefile,delimiter='\t', unpack = True,skiprows=2, usecols = (0,3))
    a[1]-=a[1,-1]
    anthracenefluorescence=RamanSpectrum(pandas.Series(a[1],a[0]))
    
     ###Normalizing to value of anthracene at 420 nm The area for the anthracene fluorescence is related to this value by 78.203    
  #  anthracenefluorescencearea = anthracenefluorescence[420]*78.2032212661
  #  anthracenefluorescencearea = anthracenefluorescence[440]*292.86
    anthracenefluorescencearea = anthracenefluorescence[470]*1257
    print 'anthracene fluorescence area=', '%.2E' % anthracenefluorescencearea
    print anthracenefluorescence.calc_area((355,550))/anthracenefluorescence[470], 'ratio of total anthracene fluorescence area to value at 470' 
    
    oneminusTdot = 1-10**(-absvalues)   ##### gives the fraction of photons absorbed by dots
    
    oneminusT_anthracene350 =1-10**(-anthraceneabsorbance350)
    print 'anthracene absorbance at 350 nm:', anthraceneabsorbance350,'. Fraction photons absorbed:', oneminusT_anthracene350
    print 'dot absorbance at 350 nm:', absvalues, '. Fraction photons absorbed:', oneminusTdot

    
    a = loadtxt(fluorescencefile,delimiter='\t', unpack = True,skiprows=2, usecols = (0,3))
    hi = RamanSpectrum(pandas.Series(a[1],a[0]))
    if subtractfluorfile !=None:
        b = loadtxt(subtractfluorfile,delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
        fluorbackground = RamanSpectrum(pandas.Series(b[1],b[0]))
        hi[:]-=fluorbackground[:]
        fluorbackground.plot(ax=fluorplot)
    hi.plot(ax=fluorplot)
    
    hi[:]-=min(hi[400:500])
    hi[:]*=0.27/(1+0.00145*158)*oneminusT_anthracene350*nliq**2/nE**2 /anthracenefluorescencearea/ oneminusTdot 
    
    
    
    
    dotfluorescencearea = hi.calc_area(fluorescencerange,fill=False)
    
    
    ## quantum yield of dots using 0.27 as QY for anthracene with o2 quenching corrrection
    print  'fluorescence (bande edg) yield of dot', dotfluorescencearea
    
    if UVVisplot is not None:
        if _plot_standard:
            anthracene.plot(ax=UVVisplot)#plot(a[0],anthracene)
        dot.plot(ax=UVVisplot,label=label)
    if fluorplot is not None:
        hi.plot(ax = fluorplot,label=label)
        if _plot_standard:
            anthracenefluorescence.plot(ax = fluorplot,label=label)
        
        
    return  dotfluorescencearea
    
def indivCdSeQY(UVVisfile, UVViscolumn, rhodaminecolumn, fluorescencefile, rhodaminefluorescencefile,
                excitationwavelength = None,standardfluorescencerange=None, baselineabsorbanceat = None,
            UVVisplot=None, fluorplot=None,
            fluorescencerange = (500,600),
            nliq=1.333,
            day=0, label=None,color = 'k'):
    print '-------------------------------------'
    print 'calculating fluorescence yield for', label,'file', fluorescencefile
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    if len(UVViscolumn)==1:
        numuvviscolumn = alphabet.find(UVViscolumn)
    elif len(UVViscolumn)==2:
        #print 'longer',alphabet.find(UVViscolumn[0]),alphabet.find(UVViscolumn[1])
        numuvviscolumn = (alphabet.find(UVViscolumn[0])+1)*26+alphabet.find(UVViscolumn[1])
    
    
    a = loadtxt(UVVisfile, delimiter = ',', unpack = True, skiprows = 1,usecols=(0,numuvviscolumn,alphabet.find(rhodaminecolumn)))
    
    a[1:]-=transpose([a[1:,0]])
    
   
    rhodamine=RamanSpectrum(pandas.Series(a[2][::-1],a[0][::-1]))
    dot=RamanSpectrum(pandas.Series(a[1][::-1],a[0][::-1]))
    
    if baselineabsorbanceat != None:
        dot-=dot[baselineabsorbanceat]
    
    rhodamine[:]-=rhodamine[700]
    
   
    rhodamineabsorbance350= rhodamine[excitationwavelength]#(rhodamine[374]-anthracene[389])*0.6735# 
    absvalues = dot[excitationwavelength]
    
    nE = 1.359
    nQ = 1.44
    nW = 1.333  ## refractive index water    
 
    
    
    a = loadtxt(rhodaminefluorescencefile,delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
    a[1]-=a[1,-1]
    standardfluorescence=RamanSpectrum(pandas.Series(a[1],a[0]))
    
     ###Normalizing to area of rhodamine B   
    #pdb.set_trace()
    standardfluorescencearea = standardfluorescence.calc_area(standardfluorescencerange)
    print 'standard fluorescence area=', '%.2E' % standardfluorescencearea
    
    oneminusTdot = 1-10**(-absvalues)   ##### gives the fraction of photons absorbed by dots
    
    oneminusT_rhodamine350 =1-10**(-rhodamineabsorbance350)  ##### gives the fraction of photons absorbed by standard
    print 'rhodamine absorbance at',excitationwavelength,'nm:', rhodamineabsorbance350,'. Fraction photons absorbed:', oneminusT_rhodamine350
    print 'dot absorbance at ',excitationwavelength,' nm:', absvalues, '. Fraction photons absorbed:', oneminusTdot

    
    a = loadtxt(fluorescencefile,delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
    hi = RamanSpectrum(pandas.Series(a[1],a[0]))
    hi[:]-=min(hi)
    hi[:]*=0.65*oneminusT_rhodamine350*nliq**2/nE**2 /standardfluorescencearea/ oneminusTdot 
    
    
    dotfluorescencearea = hi.calc_area(fluorescencerange,fill=False)
    
    
   
    print  'fluorescence (bande edge) yield of dot', dotfluorescencearea
    
    if UVVisplot is not None:
        #rhodamine.plot(ax=UVVisplot))
        dot.plot(ax=UVVisplot,label=label)
    if fluorplot is not None:
        hi.plot(ax = fluorplot,label=label)
       # standardfluorescence.plot(ax = fluorplot,label=label)
        
        
    return  dotfluorescencearea