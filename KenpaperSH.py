# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 14:32:29 2015

@author: chris
"""
from scipy.optimize import curve_fit
from ramanTools.RamanSpectrum import *
from numpy import *
from matplotlib.pyplot import *
from copy import deepcopy
import copy

from scipy.optimize import minimize
from matplotlib import gridspec

    

def SH():
    
    
    
    os.chdir('/home/chris/Documents/DataWeiss/150728')
    cdmbt=copy.deepcopy(CdMethylTPRef)
    mbt=copy.deepcopy(MethylTPRef)
    cdmbt.autobaseline((193,4000),order = 0)
    mbt.autobaseline((193,4000),order = 0)

    mbt[:]/=95
    cdmbt[:]/=5
    cdmbt.to_csv('/home/chris/Dropbox/Ken/CdMBT2.csv')
   
    
    
    fig1 = figure(figsize=(6, 12)) 
    

    
   
    
    A = RamanSpectrum('filesA.txt')  ##450 eq
    B = RamanSpectrum('filesB.txt')  #200 eq MBT
    C = RamanSpectrum('filesC.txt') #100 eq MBT
    D = RamanSpectrum('filesD.txt') # 80 eq MBT
    E = RamanSpectrum('filesE.txt') # 50 eq MBT
    F = RamanSpectrum('filesF.txt')  #25 eq MBT
    G = RamanSpectrum('/home/chris/Documents/DataWeiss/150408/150408_03.txt')
    
    G.autobaseline((2500,2700,3100,3200),order = 2, specialoption='points')
    
    for z in [A,B,C,D,E,F]:
        #z = SPIDcorrect633(z)
        z.autobaseline((200,361),order = 1,join='start')
        z.autobaseline((361,394),order = 2,join='start')
        z.autobaseline((394,647),order = 2,join='start')
        z.autobaseline((647,682),order = 0,join='start')
        z.autobaseline((682,923),order = 0,join='start')
        z.autobaseline((923,955),order = 0,join='start')
        z.autobaseline((955,1187),order = 0,join='start')
        z.autobaseline((1187,1214),order = 0,join='start')
        z.autobaseline((1214,1437),order = 0,join='start')
        z.autobaseline((1437,1462),order = 0,join='start')
        z.autobaseline((1462,1675),order = 2,join='start')
        z.autobaseline((1675,1701),order = 0,join='start')
        z.autobaseline((1701,1900), order =2,join='start')
        z.autobaseline((1900,2400), order =5,join='start')
        z.autobaseline((2400,3200),order = 0,join='start')
        z.autobaseline((3200,3600),order = 0,join='start')
       # z.autobaseline((981,1013,1098,1141,1251,1491),order =3,specialoption='points', join='start')
        
        z[:]/=50
        z.smooth()
        z-=z[2402]
        
    for z in [G]:

        z.autobaseline((2400,3200),order = 0)
        z[:]/=50
        z.smooth()
        z-=z[2402]
        
        
    mbt[:]+=700 
    mbt.set_name('mmmmmm')
    A[:]*=2
    A[:]+=600
    B[:]+=500
    C[:]+=400
    D[:]+=300
    E[:]+=200
    F[:]+=100
    G[:]*=10
    
    for z in [mbt,A,B,C,D,E,F,G]:
        z.to_csv('/home/chris/Dropbox/Ken/SHregion/'+z.name[-5]+'.csv')
        
 
    
    mbt.plot()
    A.plot()
    B.plot()
    C.plot()
    D.plot()
    E.plot()
    F.plot()
    G.plot()
    
    fs = 14
    anx = 2605
    annotate('solid MBT',(anx,740), fontsize=fs)
    annotate('450 eq MBT',(anx,620), fontsize=fs)
    annotate('200 eq MBT',(anx,520), fontsize=fs)
    annotate('100 eq MBT',(anx,420), fontsize=fs)
    annotate('80 eq MBT',(anx,320), fontsize=fs)
    annotate('50 eq MBT',(anx,220), fontsize=fs)
    annotate('25 eq MBT',(anx,120), fontsize=fs)
    annotate('0 eq',(anx,10), fontsize=fs)    
    
    xlim(2500,3030)
    ylim(-20,1300)
    return 0


    