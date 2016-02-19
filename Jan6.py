# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 12:39:39 2016

@author: chris
"""

from ramanTools.RamanSpectrum import *

def Jan6():  ### Raman spectra of PPA exchanged dots at different points in exchange
    clf()
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_03.txt')
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_05.txt')
    c = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_06.txt')
    d = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_07.txt')
    e = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_09.txt')
    f = RamanSpectrum('/home/chris/Dropbox/DataWeiss/160106/160106_10.txt')
    
    a = removespikes(a)
    a.autobaseline((200,1800), order = 5)
    a.autobaseline((1800,2100,2700,3200,3600), specialoption='points',order = 3,join='start')
    
    b.autobaseline((200,1800), order = 5)
    b.autobaseline((1800,2100,2700,3200,3600), specialoption='points',order = 3,join='start')
    b[:]+=10000
    c.autobaseline((200,1800), order = 5)
    c.autobaseline((1800,2100,2700,3200,3600), specialoption='points',order = 3,join='start')
    c[:]+=20000
    d.autobaseline((200,1800), order = 5)
    d.autobaseline((1800,2100,2700,3200,3600), specialoption='points',order = 3,join='start')
    d[:]+=30000
    e.autobaseline((200,1800), order = 5)
    e.autobaseline((1800,2100,2700,3200,3600), specialoption='points',order = 3,join='start')
    e[:]+=40000
    f.autobaseline((200,1800), order = 5)
    f.autobaseline((1800,2100,2700,3200,3600), specialoption='points',order = 3,join='start')
    f[:]+=50000
    
    a.plot(color = 'k')
    b.plot(color = 'k')
    c.plot(color = 'k')
    d.plot(color = 'k')
    e.plot(color = 'k')
    f.plot(color = 'k')
    
    OPAdots =  RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt')  ## dots
    OPAdots= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt')  ## dots
    OPAdots.autobaseline((911,1196,1385,1515,1800),specialoption='points',order=7,join='start')
    OPAdots.autobaseline((280,600,690,826,861,911),specialoption='points', order = 5,join='end')
    
    OPAdots.autobaseline((281,630), order = 4, join = 'end')    
    OPAdots.autobaseline((1492,1800), order = 5, join = 'start')
    OPAdots[:]*=3
    OPAdots[:]+=35000
   # OPAdots.plot()
    DMF = RamanSpectrum('/home/chris/Dropbox/DataWeiss/151215/151215_04.txt')
    DMF.autobaseline((300,1800), order =3)
    DMF[:]+=21000
    #DMF.plot()
    #(CdOPARef*3+35000).plot(color = 'b')
    return 0