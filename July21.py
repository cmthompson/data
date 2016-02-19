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

from scipy.optimize import minimize
from matplotlib import gridspec
def fluorescenceJuly21():
     a = loadtxt('/home/chris/Documents/DataWeiss/150721/150721_fluorescence/Stoich.csv',skiprows=1, delimiter = ',', unpack = True)
     b = loadtxt('/home/chris/Documents/DataWeiss/150721/150721_fluorescence/enriched.csv',skiprows = 1, delimiter = ',', unpack = True)
     plot(a[0],a[3])
     plot(b[0],b[3])
     def function(x,A1,w1,G1,m,b): return m*x/1000+b + A1*exp(-(x-w1)**2/G1)
     
     r = curve_fit(function, a[0],a[3],[10000,525,30,0,0])
     print r[0]
     plot(a[0],function(a[0],*r[0]))
     
     r = curve_fit(function, b[0],b[3],[10000,550,30,0,0])
     print r[0]
     plot(b[0],function(b[0],*r[0]))
        
     return 0
     
def July21UVVis():
    figure()
    ax1 = subplot(121)
    ax2 = subplot(122)
    ax1.annotate('Cd-enriched', (0.2,0.9), xycoords= 'axes fraction',fontsize = 24)
    
    a = loadtxt('/home/chris/Documents/DataWeiss/150721/150721uvvis/CdSe dots after MBT.csv',delimiter = ',', skiprows =1, unpack = True)
    x1 = argmin(abs(a[0]-535))
    x2 = argmin(abs(a[0]-545))
    xfoot = argmin(abs(a[0]-620))
    for i in a[1:7]:
        i = SGsmooth(a[0],i)

        i-=i[xfoot]
        ymax= max(i[x2:x1])
        
        i/=ymax
        ax1.plot(a[0],i)
    a = loadtxt('/home/chris/Documents/DataWeiss/150721/150721uvvis/CdSe dots before MBT.csv',delimiter = ',', skiprows =1, unpack = True)
    
    x1 = argmin(abs(a[0]-535))
    x2 = argmin(abs(a[0]-545))
    xfoot = argmin(abs(a[0]-620))
    
    
    
    a[1] = SGsmooth(a[0],i)

    a[1]-=a[1,xfoot]
    ymax= max(a[1,x2:x1])
        
    a[1]/=ymax
    ax1.plot(a[0],a[1])

    legend(['450eq','200eq','100eq','80eq','50eq','25eq','0eq'])
    
    
    ax2.annotate('Stoichiometric', (0.2,0.9), xycoords = 'axes fraction',fontsize = 24)
   
    a = loadtxt('/home/chris/Documents/DataWeiss/150721/150721uvvis/CdSe dots after MBT.csv',delimiter = ',', skiprows =1, unpack = True)
    x1 = argmin(abs(a[0]-505))
    x2 = argmin(abs(a[0]-515))
    xfoot = argmin(abs(a[0]-560))
    for i in a[array([7,8,10,9,11,12])]:
        i = SGsmooth(a[0],i)

        i-=i[xfoot]
        ymax= max(i[x2:x1])
        
        i/=ymax
        ax2.plot(a[0],i)
    a = loadtxt('/home/chris/Documents/DataWeiss/150721/150721uvvis/CdSe dots before MBT.csv',delimiter = ',', skiprows =1, unpack = True)
    
    x1 = argmin(abs(a[0]-505))
    x2 = argmin(abs(a[0]-515))
    xfoot = argmin(abs(a[0]-560))
    
    
    
    a[1] = SGsmooth(a[0],i)

    a[1]-=a[1,xfoot]
    ymax= max(a[1,x2:x1])
        
    a[1]/=ymax
    ax2.plot(a[0],a[1])

    ax2.legend(['450eq','200eq','100eq','80eq','50eq','10eq','0eq'])    
    
    
    
   
    return 0
    
def fluorescenceJuly23CdS():
     a = loadtxt('/home/chris/Documents/DataWeiss/150723/CdS-from7-22.csv',skiprows=1, delimiter = ',', unpack = True)
     plot(a[0],a[3])
   
     def function(x,A1,w1,G1,m,b): return m*x/1000+b + A1*exp(-(x-w1)**2/G1)
     
     r = curve_fit(function, a[0][20:70],a[3][20:70],[10000000,426,30,0,0])
     print r[0]
     plot(a[0][20:70],function(a[0][20:70],*r[0]))
     
     G = r[0][2]
     print 'FWHM:', 2*sqrt(-G*log(0.5))
        
     return 0
     

        
def July27():
    fig = figure(figsize=(12, 6)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
    ax1 = subplot(gs[0])
    
    A = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesA.txt')  ##450 eq
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesB.txt')  #200 eq MBT
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesC.txt') #100 eq MBT
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesD.txt') # 80 eq MBT
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesE.txt') # 50 eq MBT
    #F = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesF.txt')  #25 eq MBT
    G = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/Controlenriched.txt')  # control 0 eq
    
    X = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesX.txt')  # stoich 80 eq
    Y = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesY.txt')  # control 50 eq
    Z  = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/filesZ.txt')  # control 25 eq
    
  
    for z in [A,B,C,D,E,G,X,Y,Z]:
        #z = SPIDcorrect633(z)
        z.autobaseline((358,392),order = 0,join='start')
        z.autobaseline((392,647),order = 2,join='start')
        z.autobaseline((647,682),order = 0,join='start')
        z.autobaseline((682,923),order = 0,join='start')
        z.autobaseline((923,955),order = 0,join='start')
        z.autobaseline((955,1187),order = 0,join='start')
        z.autobaseline((1187,1214),order = 0,join='start')
        z.autobaseline((1214,1437),order = 0,join='start')
        z.autobaseline((1437,1462),order = 0,join='start')
        z.autobaseline((1462,1675),order = 0,join='start')
        z.autobaseline((1675,1701),order = 0,join='start')
        z.autobaseline((1701,1900), order = 0,join='start')
        z.autobaseline((1013,1105,1166,1248,1394,1501,1670),order =6,specialoption='points', join='start')
        z[:]-=z[1667]
        z.smooth()
        
    A[:]+=10000
    B[:]+=8000
    C[:]+=6000
    D[:]+=4000
    E[:]+=2000
    
    Z[:]+=4000
    Y[:]+=6000
    X[:]+=8000
    
     
    
    
    A.plot(color = 'k',linewidth=2)
    B.plot(color = 'k',linewidth=2)
    C.plot(color = 'k',linewidth=2)
    D.plot(color = 'k',linewidth=2)
    E.plot(color = 'k',linewidth=2)
    G.plot(color = 'k',linewidth=2)
    X.plot()
    Y.plot()
    Z.plot()
    
    

    ax1.set_xlim(500,1650)
    ax1.set_ylim(-100,12000)
    
    ax2 = subplot(gs[1])
    
    A.plot(color = 'k',linewidth=2)
    B.plot(color = 'k',linewidth=2)
    C.plot(color = 'k',linewidth=2)
    D.plot(color = 'k',linewidth=2)
    E.plot(color = 'k',linewidth=2)
    G.plot(color = 'k',linewidth=2)
    X.plot()
    Y.plot()
    Z.plot()
    xannotate = 1005
    ax2.set_xlim(1000,1110)
    ax2.set_ylim(-100,12000)
    #ax2.annotate('*', (1079,9500))
    #ax2.annotate('*', (1086,9500))
    ax2.annotate('0eq', (xannotate,240))
    ax2.annotate('50eq', (xannotate,2200))
    ax2.annotate('80eq', (xannotate,4400))
    ax2.annotate('100eq', (xannotate,6400))
    ax2.annotate('200eq', (xannotate,8400))
    ax2.annotate('450eq', (xannotate,10400))

    
    
    rlist = list()
    for z in [A,B,C,D,E]:
        r = fitspectrum(z, (1013,1104),'xGaussian', [500,500,1000,1000,1000,
                                                        1032.25,1050,  1062.17,1078,1086,
                                                        15,50,15,15,15,0,z[1160]])
#        plot(r.x,r.y,color = 'b', linewidth = 2)
        for i in r.peaks:
            plot(r.x,i)
#    g0a = 13.8   #### from july10data average of four spectra widhts
#    g0b = 21.3  ### from July10data       
#    Aratio = 0.875
#
#        
#    
#    def func(x,A0, A1,A2,A3,A4,A5,A6,A7,w1,w2,w3,w4,w5,w6,w7,G1,G2,G3,G4,G5,G6,G7,m,b): 
#        return m*x/1000+b + A0*Aratio*exp(-(x-1078.1)**2/g0a)+A0*exp(-(x-1085.5)**2/g0b)+A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2)  +A3*exp(-(x-w3)**2/G3)  + A4*exp(-(x-w4)**2/G4)  + A5*exp(-(x-w5)**2/G5)  + A6*exp(-(x-w6)**2/G6)+ A7*exp(-(x-w7)**2/G7)
#    
#    dotbound_polymerbound_ratio=list()    
#    
#    
#     plot(r.x,boundtodotpeaks)
#        
#        plot(r.x,r.y,color = 'k', linewidth = 2)
#        
#        peaks = list()
#        areas = list()
#        numpeaks = 6
#        for i in range(numpeaks):
#            A = r.params[0][1:][i]
#            x0 = r.params[0][1:][i+numpeaks]
#            G = r.params[0][1:][i+2*numpeaks]
#            x= r.x
#                
#            
#            y= A*exp(-(x-x0)**2/G)+r.params[0][1:][-2]/1000*x+r.params[0][1:][-1]
#            ar1 = A*numpy.sqrt(pi*G)
#            areas.append(ar1)
#            plot(r.x,y,'k')
#        
#        fill_between(r.x,boundtodotpeaks,m*r.x/1000+b,color = 'r')
#        ratioofdotboundareastopolymerboundarea = 100#r.params[0][0]*numpy.sqrt(pi)*(Aratio*numpy.sqrt(13.8)+numpy.sqrt(21.3))/areas[3]
#        dotbound_polymerbound_ratio.append(ratioofdotboundareastopolymerboundarea)
#    print dotbound_polymerbound_ratio
#    figure()
#    plot([1000,500,100,50],dotbound_polymerbound_ratio)
    return 0 
  

def July27UVVis():
    a = loadtxt('/home/chris/Dropbox/DataWeiss/150727/150727_UVVisCdSe after MBT.csv', delimiter = ',', skiprows = 1, unpack = True)
   
    x530 = argmin(abs(a[0]-530))
    x550=argmin(abs(a[0]-550))
    subplot(121)
    peaks=list()
    for z in a[8:]:
        z = SGsmooth(a[0], z)
        z-=min(z)
        z/=max(z[x550:x530])
        r = polyfit(a[0,x550:x530],z[x550:x530],3)
        xs = arange(550,530,-0.1)
        peaks= append(peaks, xs[argmax(polyeval(r,xs))])
        plot(a[0],z)
    legend(['0','25','50','80','100','200','450'])  
    subplot(122)
    plot([0,25,50,80,100,200,450], peaks)
   
    
    
def July28():
    """MBT series on cdse enriched dots with and without a washing step"""

    fig1 = figure(figsize=(12, 6)) 
    
    
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
    ax1 = subplot(gs[0])
    
    A = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/filesA.txt')  ##450 eq
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/filesB.txt')  #200 eq MBT
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/filesC.txt') #100 eq MBT
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/filesD.txt') # 80 eq MBT
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/filesE.txt') # 50 eq MBT
    F = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/filesF.txt')  #25 eq MBT
    G = RamanSpectrum('/home/chris/Documents/DataWeiss/150727/Controlenriched.txt')  # control 0 eq
    
    
    
    
    Aw = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/Awashed.txt')  ##450 eq with wash
    Bw = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/Bwashed.txt')  #200 eq MBT  with wash
    Cw = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/Cwashed.txt') #100 eq MBT  with wash
    Dw = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/Dwashed.txt') # 80 eq MBT  with wash
    Ew = RamanSpectrum('/home/chris/Documents/DataWeiss/150728/Ewashed.txt') # 50 eq MBT  with wash
    
  
    for z in [A,B,C,D,E,F,G,Aw,Bw,Cw,Dw,Ew]:
        #z = SPIDcorrect633(z)
        z.autobaseline((358,392),order = 0,join='start')
        z.autobaseline((392,647),order = 2,join='start')
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
        z.autobaseline((1013,1105,1166,1248,1394,1501,1670,1800),order =6,specialoption='points', join='start')
        
        z[:]-=z[1667]
        z.smooth()
        
    # make plot 1 with unwashed samples. 
    
    A[:]+=12000    
    B[:]+=10000
    C[:]+=8000
    D[:]+=6000
    E[:]+=4000
    F[:]+=2000
    
    A.plot(color = 'k',linewidth=2)
    B.plot(color = 'k',linewidth=2)
    C.plot(color = 'k',linewidth=2)
    D.plot(color = 'k',linewidth=2)
    E.plot(color = 'k',linewidth=2)
    F.plot(color = 'k',linewidth=2)
    G.plot(color = 'k',linewidth=2)
    
    

    ax1.set_xlim(500,1650)
    ax1.set_ylim(-100,14000)
    
    ax2 = subplot(gs[1])
    
    A.plot(color = 'k',ax=ax2, linewidth=2)
    B.plot(color = 'k',ax=ax2,linewidth=2)
    C.plot(color = 'k',ax=ax2,linewidth=2)
    D.plot(color = 'k',ax=ax2,linewidth=2)
    E.plot(color = 'k',ax=ax2,linewidth=2)
    F.plot(color = 'k',ax=ax2,linewidth=2)
    G.plot(color = 'k',ax=ax2,linewidth=2)
    
    xannotate = 1005
    ax2.set_xlim(1000,1110)
    ax2.set_ylim(-100,12000)
    ax2.annotate('0eq', (xannotate,240))
    ax2.annotate('25eq', (xannotate,2200))
    ax2.annotate('50eq', (xannotate,4400))
    ax2.annotate('80eq', (xannotate,6400))
    ax2.annotate('100eq', (xannotate,8400))
    ax2.annotate('200eq', (xannotate,10400))
    ax2.annotate('450eq', (xannotate,12400))


    
   # def baseline(x,k,m,b):return k*G[1013:1070]+m/1000*x+b 
        
    def baseline(x,k,A1,A2, w1,w2,G1,G2,m,b):return k*G[1013:1100] + A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2) +m*x/1000+b 
   
    rlist = list()
    guess = [  9.04030796e-01 ,  4.48825353e+02,  6.73947343e+02  , 1.07795605e+03,    1.08599115e+03 ,  1.03549466e+01 ,  3.10430178e+01, -1.69440682e+02]
    for z in [A,B,C,D,E]:
#        r = fitspectrum(z, (1013,1104),'xGaussian', [500,500,1000,1000,1000,
#                                                        1032.25,1050,  1062.17,1078,1086,
#                                                        15,50,15,15,15,0,z[1160]])
    
         r = fitspectrum(z, (1013,1100),'Custom', guess+[z[1160]],function = baseline)
         print r.params[0]
         (k,A1,A2, w1,w2,G1,G2,m,b)=r.params[0]
         plot(r.x,k*G[1013:1100]+ m*r.x/1000+b )
         plot(r.x,A1*exp(-(r.x-w1)**2/G1) +m*r.x/1000+b )
         plot(r.x,A2*exp(-(r.x-w2)**2/G2) +m*r.x/1000+b )
         plot(r.x,r.y)

#
#        for i in r.peaks:
#            plot(r.x,i)

    ########### make plot 2 with washed samples 

    
    fig = figure(figsize=(12, 6))
    
    
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
    ax1 = subplot(gs[0])
    ax2 = subplot(gs[1])
    ax1.annotate('Washed samples', (100,13000),fontsize = 24)
    
    Aw[:]+=12000    
    Bw[:]+=10000
    Cw[:]+=8000
    Dw[:]+=6000
    Ew[:]+=4000
   # Fw[:]+=2000
    
    Aw.plot(color = 'k',ax=ax1,linewidth=2)
    Bw.plot(color = 'k',ax=ax1,linewidth=2)
    Cw.plot(color = 'k',ax=ax1,linewidth=2)
    Dw.plot(color = 'k',ax=ax1,linewidth=2)
    Ew.plot(color = 'k',ax=ax1,linewidth=2)
    #Fw.plot(color = 'k',linewidth=2)
    #Gw.plot(color = 'k',linewidth=2)
    
    

    ax1.set_xlim(500,1650)
    ax1.set_ylim(-100,14000)
    
    
    
    Aw.plot(color = 'k',ax=ax2,linewidth=2)
    Bw.plot(color = 'k',ax=ax2,linewidth=2)
    Cw.plot(color = 'k',ax=ax2,linewidth=2)
    Dw.plot(color = 'k',ax=ax2,linewidth=2)
    Ew.plot(color = 'k',ax=ax2,linewidth=2)
    #Fw.plot(color = 'k,ax=ax2',linewidth=2)
    #Gw.plot(color = 'k,ax=ax2',linewidth=2)
    xannotate = 1005
    ax2.set_xlim(1000,1110)
    ax2.set_ylim(-100,12000)
    ax2.annotate('0eq', (xannotate,240))
    ax2.annotate('25eq', (xannotate,2200))
    ax2.annotate('50eq', (xannotate,4400))
    ax2.annotate('80eq', (xannotate,6400))
    ax2.annotate('100eq', (xannotate,8400))
    ax2.annotate('200eq', (xannotate,10400))
    ax2.annotate('450eq', (xannotate,12400))

        
    #def baseline(x,k,A1,A2, w1,w2,G1,G2,m,b):return k*G[1013:1100] + A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2) +m*x/1000+b 
   
    rlist = list()
    guess = [  9.04030796e-01 ,  4.48825353e+02,  6.73947343e+02  , 1.07795605e+03,    1.08599115e+03 ,  1.03549466e+01 ,  3.10430178e+01, -1.69440682e+02]
    for z in [Aw,Bw,Cw,Dw,Ew]:
        pass
#        r = fitspectrum(z, (1013,1104),'xGaussian', [500,500,1000,1000,1000,
#                                                        1032.25,1050,  1062.17,1078,1086,
#                                                        15,50,15,15,15,0,z[1160]])
    
#         r = fitspectrum(z, (1013,1100),'Custom', guess+[z[1160]],function = baseline)
#         print r.params[0]
#         (k,A1,A2, w1,w2,G1,G2,m,b)=r.params[0]
#         plot(r.x,k*G[1013:1100]+ m*r.x/1000+b )
#         plot(r.x,A1*exp(-(r.x-w1)**2/G1) +m*r.x/1000+b )
#         plot(r.x,A2*exp(-(r.x-w2)**2/G2) +m*r.x/1000+b )
#         plot(r.x,r.y)

#
#        for i in r.peaks:
#            plot(r.x,i)
    return 0
         



def July30():  
    """MBT series on cdse enriched dots with and without a washing step"""
    os.chdir('/home/chris/Documents/DataWeiss/150730')
    fig1 = figure(figsize=(12, 6)) 
    
    
#    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
#    ax1 = subplot(gs[0])
    
    A = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/A.txt')  ##450 eq
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/B.txt')  #200 eq MBT
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/C.txt') #100 eq MBT
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/D.txt') # 80 eq MBT
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/E.txt') # 50 eq MBT
    F = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/F.txt')  #25 eq MBT
    G = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/G.txt')  # control 0 eq
    
    
    
    Aw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Awash.txt')  ##450 eq with wash
    Bw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Bwash.txt')  #200 eq MBT  with wash
    Bw[:]/=4
    Cw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Cwash.txt') #100 eq MBT  with wash
    Dw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Dwash.txt') # 80 eq MBT  with wash
    Ew = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Ewash.txt') # 50 eq MBT  with wash
    Fw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Fwash.txt') # 25 eq MBT  with wash
    Gw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Gwash.txt') # 0 eq MBT  with wash
    
    U = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/U.txt')  ##450 eq
    V = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/V.txt')  #200 eq MBT
    W = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/W.txt') #100 eq MBT
    X = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/X.txt') # 80 eq MBT
    Y = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Y.txt') # 50 eq MBT
    Z = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Z.txt')  #25 eq MBT
    ZZ = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/ZZ.txt')  # control 0 eq
    
    
    
    Uw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Uwash.txt')  ##450 eq with wash
    Vw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Vwash.txt')  #200 eq MBT  with wash
    Ww = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Wwash.txt')  #100 eq MBT  with wash
    Xw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Xwash.txt') #80 eq MBT  with wash
    Yw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Ywash.txt') # 50 eq MBT  with wash
    Zw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/Zwash.txt') # 25 eq MBT  with wash
    ZZw = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/ZZwash.txt') # 0 eq MBT  with wash


    
  
    for z in [A,B,C,D,E,F,G,Aw,Bw,Cw,Dw,Ew,Fw,Gw,U,V,W,X,Y,Z,ZZ,Uw,Vw,Ww,Xw,Yw,Zw,ZZw]:

        z.autobaseline((878,1173),order = 2)
        z.autobaseline((1142,1173),order = 0,join='start')
        z.autobaseline((1173,1300),order = 0,join='start')
        z.autobaseline((878,993,1054,1156,1201),order =4,specialoption='points', join='start')
#        
        z[:]-=z[1120]
        z.smooth()


    
    ax2 = subplot(141)#subplot(gs[1])
    ax2.annotate('Unwashed enriched', (0.05,0.8), xycoords = 'axes fraction',fontsize = 24)
    
    A.plot(color = 'r',ax=ax2, linewidth=2)
    B.plot(color = 'm',ax=ax2,linewidth=2)
    C.plot(color = 'y',ax=ax2,linewidth=2)
    D.plot(color = 'g',ax=ax2,linewidth=2)
    E.plot(color = 'b',ax=ax2,linewidth=2)
    F.plot(color = 'c',ax=ax2,linewidth=2)
    G.plot(color = 'k',ax=ax2,linewidth=2)
    quickoffset(ax2,rnge = (1000,1140))
    
    xannotate = 1005
    ax2.set_xlim(1000,1110)
    ax2.set_ylim(-100,12000)
    ax2.annotate('0eq', (xannotate,240))
    ax2.annotate('25eq', (xannotate,2200))
    ax2.annotate('50eq', (xannotate,4400))
    ax2.annotate('80eq', (xannotate,6400))
    ax2.annotate('100eq', (xannotate,8400))
    ax2.annotate('200eq', (xannotate,10400))
    ax2.annotate('450eq', (xannotate,12400))
            

    ax2 = subplot(142)
    ax2.annotate('Washed enriched', (0.05,0.8),xycoords = 'axes fraction',fontsize = 24)

    
    Aw.plot(color = 'r',ax=ax2,linewidth=2)
    Bw.plot(color = 'm',ax=ax2,linewidth=2)
    Cw.plot(color = 'y',ax=ax2,linewidth=2)
    Dw.plot(color = 'g',ax=ax2,linewidth=2)
    Ew.plot(color = 'b',ax=ax2,linewidth=2)
    
    Fw.plot(color = 'c',ax=ax2,linewidth=2)
    Gw.plot(color = 'k',ax=ax2,linewidth=2)
    quickoffset(ax2,rnge = (1000,1140))    
    
    xannotate = 1005
    ax2.set_xlim(1000,1110)
    ax2.set_ylim(-100,12000)
    ax2.annotate('0eq', (xannotate,240))
    ax2.annotate('25eq', (xannotate,2200))
    ax2.annotate('50eq', (xannotate,4400))
    ax2.annotate('80eq', (xannotate,6400))
    ax2.annotate('100eq', (xannotate,8400))
    ax2.annotate('200eq', (xannotate,10400))
    ax2.annotate('450eq', (xannotate,12400))



    ax2 = subplot(143)
    U.plot(color = 'r',ax=ax2,linewidth=2)
    V.plot(color = 'm',ax=ax2,linewidth=2)
    W.plot(color = 'y',ax=ax2,linewidth=2)
    X.plot(color = 'g',ax=ax2,linewidth=2)
    Y.plot(color = 'b',ax=ax2,linewidth=2)
    Z.plot(color = 'c',ax=ax2,linewidth=2)
    ZZ.plot(color = 'k',ax=ax2,linewidth=2)
    quickoffset(ax2,rnge = (1000,1140))    
    
    xannotate = 1005
    ax2.set_xlim(1000,1110)
    ax2.set_ylim(-100,12000)
    ax2.annotate('0eq', (xannotate,240))
    ax2.annotate('25eq', (xannotate,2200))
    ax2.annotate('50eq', (xannotate,4400))
    ax2.annotate('80eq', (xannotate,6400))
    ax2.annotate('100eq', (xannotate,8400))
    ax2.annotate('200eq', (xannotate,10400))
    ax2.annotate('450eq', (xannotate,12400))        
    

    ax2 = subplot(144)# subplot(gs[1])
    ax2.annotate('Washed stoichiometric', (0.05,0.8),xycoords = 'axes fraction', fontsize = 24)
    


    
    Uw.plot(color = 'r',ax=ax2,linewidth=2)
    Vw.plot(color = 'm',ax=ax2,linewidth=2)
    Ww.plot(color = 'y',ax=ax2,linewidth=2)
    Xw.plot(color = 'g',ax=ax2,linewidth=2)
    Yw.plot(color = 'b',ax=ax2,linewidth=2)
    Zw.plot(color = 'c',ax=ax2,linewidth=2)
    ZZw.plot(color = 'k',ax=ax2,linewidth=2)
    quickoffset(ax2,rnge = (1000,1140))    
    
    xannotate = 1005
    ax2.set_xlim(1000,1110)
    ax2.set_ylim(-100,12000)
    ax2.annotate('0eq', (xannotate,240))
    ax2.annotate('25eq', (xannotate,2200))
    ax2.annotate('50eq', (xannotate,4400))
    ax2.annotate('80eq', (xannotate,6400))
    ax2.annotate('100eq', (xannotate,8400))
    ax2.annotate('200eq', (xannotate,10400))
    ax2.annotate('450eq', (xannotate,12400))
    
    
    (i.set_xlim(1040,1140) for i in gcf().get_axes()) 
            
    return 0
    
def Aug3():
    
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150803/sample U/sample washed stoich_3.txt') # 0 eq MBT  with wash
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150803/sample U/sample washed stoich_2.txt') # 0 eq MBT  with wash
    c = RamanSpectrum('/home/chris/Documents/DataWeiss/150803/sample U/sample washed stoich_4.txt') # 0 eq MBT  with wash


    
  
    for z in [a,b,c]:

        z.autobaseline((900,1180),order = 2)
        
      #  z.autobaseline((878,993,1054,1156,1201),order =4,specialoption='points', join='start')
#        
        z[:]-=z[1120]
        z.smooth(window_len = 11, window = 'SG')
        z.plot()
    return 0
    
def FinalKen():
    """data finally used for Ken's paper.  MBT treated dots for stoichiometric and Cd-enriched"""

    os.chdir('/home/chris/Documents/DataWeiss/150730')
    fig1 = figure(figsize=(12, 6)) 

    
    A = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/A.txt')  ##450 eq
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/B.txt')  #200 eq MBT
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/C.txt') #100 eq MBT
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/D.txt') # 80 eq MBT
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/E.txt') # 50 eq MBT
    F = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/F.txt')  #25 eq MBT
    G = RamanSpectrum('/home/chris/Documents/DataWeiss/150730/G.txt')  # control 0 eq

    for z in [A,B,C,D,E,F,G]:

        z.autobaseline((878,1173),order = 2)
        z.autobaseline((1142,1173),order = 0,join='start')
        z.autobaseline((1173,1300),order = 0,join='start')
        z.autobaseline((878,993,1054,1104,1156,1201),order =4,specialoption='points', join='start')
#        
        z[:]-=z[1105]
        z.smooth()


   
    ax2 = subplot(121)#subplot(gs[1])
    ax2.annotate('Cd-Enriched',  (0.95,0.95), xycoords = 'axes fraction',fontsize = 18,horizontalalignment = 'right')
    
    
    A[:]+=6000
    B[:]+=4800
    C[:]+=4000
    D[:]+=3200
    E[:]+=2400 
    F[:]+=1200
    
    A.plot(color = 'r',ax=ax2, linewidth=2)
    B.plot(color = 'm',ax=ax2,linewidth=2)
    C.plot(color = 'y',ax=ax2,linewidth=2)
    D.plot(color = 'g',ax=ax2,linewidth=2)
    E.plot(color = 'b',ax=ax2,linewidth=2)
    F.plot(color = 'c',ax=ax2,linewidth=2)
    G.plot(color = 'k',ax=ax2,linewidth=2)
    
    
    g0a = 20   #### from july10data average of four spectra widhts
    g0b = 12  ### from July10data       
    Aratio = 1.2

        
    
    def func(x,A0, A1,A2,A3,A4,w1,w2,w3,w4,G1,G2,G3,G4,m,b): 
        
        return m*x/1000+b + A0*Aratio*exp(-(x-1079)**2/g0a)+A0*exp(-(x-1086)**2/g0b)+A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2)  +A3*exp(-(x-w3)**2/G3)  + A4*exp(-(x-w4)**2/G4) #+ A5*exp(-(x-w5)**2/G5)#  + A6*exp(-(x-w6)**2/G6)+ A7*exp(-(x-w7)**2/G7)
    
    dotbound_polymerbound_ratio=list()  
    
#    for z in [A]:
#        guess = [345.29413223,   262.03951404,   900, 800,1000, 140.81815835,
#                         1063.37326114,  1062.05313338, 1079.11296399,  1086.09421142, 1089,1090.78457976,
#                                 54.25792561,     5.98554751 ,  15  ,   12.71894669 ,15  , 35.70134452,
#                                     0,    z[1100]]
#        print guess[-1]
#        r = fitspectrum(z, (1050,1100),'xGaussian',guess)
# 
#        
#        for p in r.peaks:
#            plot(r.x,p,'k')
#        plot(r.x,r.y,'k')
#    for z in [B]:
#        guess = [345.29413223,   262.03951404,   404.93428644,   215.18470604 ,200, 140.81815835,
#                         1063.37326114,  1062.05313338, 1079.11296399,  1086.09421142, 1089,1090.78457976,
#                                 54.25792561,     5.98554751 ,  32.9496268  ,   12.71894669 ,15  , 35.70134452,
#                                     -10    ,z[1100]]
#        print guess[-1]
#        r = fitspectrum(z, (1050,1100),'xGaussian',guess)
#        print r.params[0][6:9]
#        
#        
#        for p in r.peaks:
#            plot(r.x,p,'k')
#        plot(r.x,r.y,'k')
    for z in [G]:
        
        r = fitspectrum(z, (1050,1100),'xGaussian', [200,700,700,700,100,
                                                          1062,1061.76,  1078.21,  1085.099,1095,
                                                      15,15,15,15,15,0,z[1100]])
        print r.params[0]
        print 'ratio of A',  r.params[0][1]/r.params[0][2]
        print 'g0,g1 = ',r.params[0][12],r.params[0][13]
        m = r.params[0][-2]
        b =r.params[0][-1]
        for p in range(len(r.peaks)):

                plot(r.x,r.peaks[p],'k')
        plot(r.x,r.y,'k')
    
    for z in [C,D,E,F]:
        
        r = fitspectrum(z, (1050,1100),'xGaussian', [200,700,700,700,100,
                                                          1062,1061.76,  1078.21,  1085.099,1095,
                                                      15,15,15,15,15,0,z[1100]])
        print r.params[0]
        print 'ratio of A',  r.params[0][1]/r.params[0][2]
        print 'g0,g1 = ',r.params[0][12],r.params[0][13]
        m = r.params[0][-2]
        b =r.params[0][-1]
        for p in range(len(r.peaks)):
            if p ==2:
                fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='g')
            elif p==3:
                fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='b')
            else:
                plot(r.x,r.peaks[p],'k')
        plot(r.x,r.y,'k')
    print 'bbbbbbbbbbbbbb'
    for z in [B]:
     
        guess = [1000,
                         345.29413223,   262.03951404,   100, 100,
                         1063.37326114,  1062.05313338,  1089,1097,
                                 54.25792561,     5.98554751 ,15  ,7, 
                                     -2000,    z[1100]+1.1*2000]
        r = fitspectrum(z, (1050,1103),'Custom', guess,function=func)
        print r.params[0]
        A0 = r.params[0][0]
        m = r.params[0][-2]
        b =r.params[0][-1]
        boundtodotpeak1 = A0*Aratio*exp(-(r.x-1079)**2/g0a)+m*r.x/1000+b
        boundtodotpeak2 = A0*exp(-(r.x-1086)**2/g0b)+m*r.x/1000+b
        
        fill_between(r.x,boundtodotpeak1,m*r.x/1000+b,color ='g')
        fill_between(r.x,boundtodotpeak2,m*r.x/1000+b,color ='b')
        plot(r.x,r.y,color = 'k', linewidth = 2)
        
        peaks = list()
        areas = list()
        numpeaks = (len(guess)-3)/3
        print 'num peaks detected', numpeaks
        for i in range(numpeaks):
            
            A0 = r.params[0][1:][i]
            x0 = r.params[0][1:][i+numpeaks]
            G0 = r.params[0][1:][i+2*numpeaks]
            m = r.params[0][1:][-2]
            b = r.params[0][1:][-1]
            x= r.x
                
            
            y= A0*exp(-(x-x0)**2/G0)+m*x/1000+b
           
            ar1 = A0*numpy.sqrt(pi*G0)
            areas.append(ar1)
            if i == argmin(abs(r.params[0][5:9]-1089)):
                fill_between(r.x,y,m*r.x/1000+b,color ='r')
            else:
                ax2.plot(r.x,y,'k')
    for z in [A]:
        guess = [500,
                         345.29413223,   262.03951404,   500, 10,
                         1063.37326114,  1062.05313338,  1089,1097,
                                 54.25792561,     5.98554751 ,15  ,7, 
                                     -2000,    z[1100]+1.100*2000]
        r = fitspectrum(z, (1050,1103),'Custom', guess,function=func)
        print r.params[0]
        A0 = r.params[0][0]
        m = r.params[0][-2]
        b =r.params[0][-1]
        boundtodotpeak1 = A0*Aratio*exp(-(r.x-1079)**2/g0a)+m*r.x/1000+b
        boundtodotpeak2 = A0*exp(-(r.x-1086)**2/g0b)+m*r.x/1000+b
        
        fill_between(r.x,boundtodotpeak1,m*r.x/1000+b,color ='g')
        fill_between(r.x,boundtodotpeak2,m*r.x/1000+b,color ='b')
        plot(r.x,r.y,color = 'k', linewidth = 2)
        
        peaks = list()
        areas = list()
        numpeaks = (len(guess)-3)/3
        print 'num peaks detected', numpeaks
        for i in range(numpeaks):
            A0 = r.params[0][1:][i]
            x0 = r.params[0][1:][i+numpeaks]
            G0 = r.params[0][1:][i+2*numpeaks]
            m = r.params[0][1:][-2]
            b = r.params[0][1:][-1]
            x= r.x
                
            
            y= A0*exp(-(x-x0)**2/G0)+m*x/1000+b
           
            ar1 = A0*numpy.sqrt(pi*G0)
            areas.append(ar1)
            if i == argmin(abs(r.params[0][5:9]-1089.3)):
                fill_between(r.x,y,m*r.x/1000+b,color ='r')
            else:
                ax2.plot(r.x,y,'b-')
            
#        ratioofdotboundareastopolymerboundarea = r.params[0][0]*numpy.sqrt(pi)*(Aratio*numpy.sqrt(13.8)+numpy.sqrt(21.3))/areas[3]
 #       dotbound_polymerbound_ratio.append(ratioofdotboundareastopolymerboundarea)

        #fill_between(r.x,boundtodotpeaks,m*r.x/1000+b,color ='r')



    
    xannotate = 1155
    yup = 100
    
    ax2.annotate('0eq', (xannotate,G[1160]+yup),horizontalalignment='right')
    ax2.annotate('25eq', (xannotate,F[1160]+yup),horizontalalignment='right')
    ax2.annotate('50eq', (xannotate,E[1160]+yup),horizontalalignment='right')
    ax2.annotate('80eq', (xannotate,D[1160]+yup),horizontalalignment='right')
    ax2.annotate('100eq', (xannotate,C[1160]+yup),horizontalalignment='right')
    ax2.annotate('200eq', (xannotate,B[1160]+yup),horizontalalignment='right')
    ax2.annotate('450eq', (xannotate,A[1160]+yup),horizontalalignment='right')
    ax2.set_xlim(1040,1160)
    ax2.set_ylim(-100,8000)
    
    
    
    
     
    ax2 = subplot(122)#subplot(gs[1])
    ax2.annotate('Stoichiometric', (0.95,0.95), xycoords = 'axes fraction',fontsize = 18,horizontalalignment = 'right')
    
    print 'stoichiometric fits:'
  
   
    
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_02.txt')
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_03.txt')
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_04.txt')
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_05.txt')
   
    
  
    for z in [B,C,D,E]:
        
        z.autobaseline((358,392),order = 1,join='start')
        z.autobaseline((392,647),order = 2,join='start')
        z.autobaseline((647,682),order = 1,join='start')
        z.autobaseline((682,923),order = 2,join='start')
        z.autobaseline((923,955),order = 1,join='start')
        z.autobaseline((955,1187),order = 2,join='start')
        z.autobaseline((1187,1214),order = 1,join='start')
        z.autobaseline((1214,1437),order = 2,join='start')
        z.autobaseline((1437,1462),order = 1,join='start')
        z.autobaseline((1462,1675),order = 1,join='start')
        z.autobaseline((1675,1701),order = 1,join='start')
        z.autobaseline((1701,1900), order = 1,join='start')
        z.autobaseline((1701,1900), order = 1,join='start')        
        
        z[:]-=z[1050]
        z[:]*=0.7
        z.smooth()
    

    B[:]+=6200 ##500 eq
    C[:]+=4200 ##100 eq
    D[:]+=2600 ##50eq
    E[:]*=3
    E[:]+=700
    

    
#    (ODPARef/6+8000).plot()
#    (MethylTPRef/20+8000).plot()
#    (CdMethylTPRef/2+6000).plot()

  
    B.plot(color = 'r',linewidth=2)
    C.plot(color = 'y',linewidth=2)
    D.plot(color = 'b',linewidth=2)
    E.plot(color = 'k',linewidth=2)
    
    ax2.set_xlim(1040,1160)
    ax2.set_ylim(-10,8000)
   
    ax2.annotate('0eq', (1130,1020))
    ax2.annotate('50eq', (1130,2930))
    ax2.annotate('100eq', (1130,4880))
    ax2.annotate('500eq', (1130,6880))
    
    for z in [B,C,D] :
        r = fitspectrum(z, (1047,1100),'xGaussian', [500,1500,1500,
                                                         1064.76,  1078.21,  1085.099,
                                                        50,15,15,0,z[1100]])
        print r.params[0]
       
        m = r.params[0][-2]
        b =r.params[0][-1]
        for p in range(len(r.peaks)):
            if p ==1:
                fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='g')
            elif p==2:
                fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='b')
            else:
                plot(r.x,r.peaks[p],'k')
        plot(r.x,r.y,'k')
    for z in [E] :
        r = fitspectrum(z, (1047,1100),'xGaussian', [500,1500,1500,
                                                         1064.76,  1078.21,  1085.099,
                                                        50,15,15,0,z[1100]])
        print r.params[0]
       
        m = r.params[0][-2]
        b =r.params[0][-1]
        for p in range(len(r.peaks)):

            plot(r.x,r.peaks[p],'k')
        plot(r.x,r.y,'k')
    
        
    

            
            
            
    return 0
    



    