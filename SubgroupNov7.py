# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 08:48:21 2014

@author: chris

Subgroup files for Nov7
"""
from ramanTools import RamanSpectrum as RamanTools


def Fig1():
    figure()
    a = RamanTools.RamanSpectrum('/home/chris/Documents/DataWeiss/141105/18_1.txt')
    
    
    a.autobaseline((1050,1800))
    #a._smooth(window_len=7)
    a[:]/=5
    
  
    
    
    b = RamanTools.RamanSpectrum('/home/chris/Documents/DataWeiss/141105/20_1.txt')
    b[:]-=1450
    #b._smooth(window_len=7)
    b.autobaseline((860,1800))
   
   
    
    c= RamanTools.RamanSpectrum('/home/chris/Documents/DataWeiss/141105/21_1.txt')
    c.autobaseline((860,1800))
    #c._smooth(window_len=7)
    
    
    
    b+=100
    c+=200
    
    a.plot()
    b.plot()    
    c.plot() 
    
    ylim(-100,800)
 
    
    annotate('*',(1206,400))
    annotate('*',(1448,680))
    annotate('*: TCNQ0',(1602,440))
    #annotate('*: TCNQ0',(0.9,0.99),textcoords = 'axes fraction')
    legend(['PbS+10TCNQ','PbS+30TCNQ','PbS+30TCNQ showing TCNQ0'],loc = 2)
    ylabel('Intensity (a.u.)')
    xlabel('Raman shift (cm$^{-1}$)')
    return 0
    

    
    