# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 08:38:59 2014

@author: chris
"""
import RamanTools2 as RamanTools

def Fig1():  ### show the phonon
    for filename in [#'/home/chris/Documents/DataWeiss/141113/1_10TCNQ on CdSE_1.txt',
                     #'/home/chris/Documents/DataWeiss/141113/2_10TCNQ on CdSE_1.txt',
                     '/home/chris/Documents/DataWeiss/141113/3_1.txt']:
                     #'/home/chris/Documents/DataWeiss/141113/10_1.txt']:
        a = RamanTools.RamanSpectrum(filename)
        
      
        a.plot(color = 'k')

    return 0
    