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
import weissdatavariables
    
def FinalKen():
    
    
    
    os.chdir('/home/chris/Dropbox/DataWeiss/150730')
    cdmbt=copy.copy(CdMethylTPRef)
    mbt=copy.copy(MethylTPRef)
    cdmbt.autobaseline((1000,1200),order = 0)
    mbt.autobaseline((1000,1200),order = 0)

    mbt[:]/=95
    cdmbt[:]/=5
    
    
    
    fig1 = figure(figsize=(6, 12)) 
    
    gs = gridspec.GridSpec(2, 2, height_ratios=[1,5]) 
    ax1top = subplot(gs[0]) 
    ax2top = subplot(gs[1])
    ax1 = subplot(gs[2])
    ax2 =  subplot(gs[3])
    
   
    
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



   
    #ax2 = subplot(gs[2])#subplot(gs[1])
    ax1.annotate('Cd-Enriched',  (0.95,0.95), xycoords = 'axes fraction',fontsize = 18,horizontalalignment = 'right')
    
    
    A[:]+=6000
    B[:]+=4800
    C[:]+=4000
    D[:]+=3200
    E[:]+=2400 
    F[:]+=1200
    
    A.plot(color = 'r',ax=ax1, linewidth=2)
    B.plot(color = 'm',ax=ax1,linewidth=2)
    C.plot(color = 'y',ax=ax1,linewidth=2)
    D.plot(color = 'g',ax=ax1,linewidth=2)
    E.plot(color = 'b',ax=ax1,linewidth=2)
    F.plot(color = 'c',ax=ax1,linewidth=2)
    G.plot(color = 'k',ax=ax1,linewidth=2)
    return 0
    
    g0a = 20   #### from july10data average of four spectra widhts
    g0b = 12  ### from July10data       
    Aratio = 1.2

        
    
    def func(x,A0, A1,A2,A3,A4,w1,w2,w3,w4,G1,G2,G3,G4,m,b): 
        
        return m*x/1000+b + A0*Aratio*exp(-(x-1079)**2/g0a)+A0*exp(-(x-1086)**2/g0b)+A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2)  +A3*exp(-(x-w3)**2/G3)  + A4*exp(-(x-w4)**2/G4) #+ A5*exp(-(x-w5)**2/G5)#  + A6*exp(-(x-w6)**2/G6)+ A7*exp(-(x-w7)**2/G7)
    
    dotbound_polymerbound_ratio=list()  
    

    for z in [G]:
        
        
        r = fitspectrum(z, (1050,1100),'xGaussian', [200,700,700,700,100,
                                                          1062,1061.76,  1078.21,  1085.099,1095,
                                                      15,15,15,15,15,0,z[1100]])
                                                      
        tosave = array([r.x])
        tosave = append(tosave,[r.y_in],axis=0)        
        tosave = append(tosave,[r.y],axis=0)
        tosave=append(tosave,r.peaks,axis=0)
        tosave=transpose(tosave)
        #savetxt('/home/chris/Ken/Cdenriched/'+z.name[-5]+'.csv',tosave,header='x,baselined data,total fit'+(tosave.shape[1]-3)*',peak')
        
        print r.params[0]
        print 'ratio of A',  r.params[0][1]/r.params[0][2]
        print 'g0,g1 = ',r.params[0][12],r.params[0][13]
        m = r.params[0][-2]
        b =r.params[0][-1]
        for p in range(len(r.peaks)):

                ax1.plot(r.x,r.peaks[p],'k')
        ax1.plot(r.x,r.y,'k')
    
    for z in [C,D,E,F]:
        
        r = fitspectrum(z, (1050,1100),'xGaussian', [200,700,700,700,100,
                                                          1062,1061.76,  1078.21,  1085.099,1095,
                                                      15,15,15,15,15,0,z[1100]])
        print r.params[0]
        print 'ratio of A',  r.params[0][1]/r.params[0][2]
        print 'g0,g1 = ',r.params[0][12],r.params[0][13]
        m = r.params[0][-2]
        b =r.params[0][-1]
        
        tosave = array([r.x])
        tosave = append(tosave,[r.y_in],axis=0)        
        tosave = append(tosave,[r.y],axis=0)
        tosave=append(tosave,r.peaks,axis=0)
        tosave=transpose(tosave)
        #savetxt('/home/chris/Ken/Cdenriched/'+z.name[-5]+'.csv',tosave,delimiter=',',header='x,baselined data,total fit'+(tosave.shape[1]-3)*',peak')
        for p in range(len(r.peaks)):
           #tosave = append(tosave,r.peaks[p],axis=1)
            if p ==2:
                ax1.fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='g')
            elif p==3:
                ax1.fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='b')
            else:
                ax1.plot(r.x,r.peaks[p],'k')
        ax1.plot(r.x,r.y,'k')
        
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
        
        tosave = array([r.x])
        tosave = append(tosave,[r.y_in],axis=0)        
        tosave = append(tosave,[r.y],axis=0)
        
     
        
        
        
    
        
        ax1.fill_between(r.x,boundtodotpeak1,m*r.x/1000+b,color ='g')
        ax1.fill_between(r.x,boundtodotpeak2,m*r.x/1000+b,color ='b')
        ax1.plot(r.x,r.y,color = 'k', linewidth = 2)
        
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
            tosave=append(tosave,[y],axis=0)
            if i == argmin(abs(r.params[0][5:9]-1089)):
                ax1.fill_between(r.x,y,m*r.x/1000+b,color ='r')
            else:
                ax1.plot(r.x,y,'k')
        tosave=append(tosave,[boundtodotpeak1],axis=0)
        tosave=append(tosave,[boundtodotpeak2],axis=0)
        tosave=transpose(tosave)
        #savetxt('/home/chris/Ken/Cdenriched/'+z.name[-5]+'.csv',tosave,delimiter=',',header='x,baselined data,total fit'+(tosave.shape[1]-3)*',peak')
    for z in [A]:
        guess = [500,
                         345.29413223,   262.03951404,   500, 10,
                         1063.37326114,  1062.05313338,  1089,1097,
                                 54.25792561,     5.98554751 ,15  ,7, 
                                     -2000,    z[1100]+1.100*2000]
        r = fitspectrum(z, (1050,1103),'Custom', guess,function=func)
        
        tosave = array([r.x])
        tosave = append(tosave,[r.y_in],axis=0)        
        tosave = append(tosave,[r.y],axis=0)
        
        A0 = r.params[0][0]
        m = r.params[0][-2]
        b =r.params[0][-1]
        boundtodotpeak1 = A0*Aratio*exp(-(r.x-1079)**2/g0a)+m*r.x/1000+b
        boundtodotpeak2 = A0*exp(-(r.x-1086)**2/g0b)+m*r.x/1000+b
        btp1_nobase=A0*Aratio*exp(-(r.x-1079)**2/g0a)
        btp2_nobase= A0*exp(-(r.x-1086)**2/g0b)
        
        tosave = array([r.x])
        tosave = append(tosave,[r.y_in],axis=0)        
        tosave = append(tosave,[r.y],axis=0)
        
        
        ax1.fill_between(r.x,boundtodotpeak1,m*r.x/1000+b,color ='g')
        ax1.fill_between(r.x,boundtodotpeak2,m*r.x/1000+b,color ='b')
        ax1.plot(r.x,r.y,color = 'k', linewidth = 2)
        
        
        ax1top.fill_between(r.x,btp1_nobase,0,color ='g')
        ax1top.fill_between(r.x,btp2_nobase,0,color ='b')
        
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
            tosave=append(tosave,[y],axis=0)
            if i == argmin(abs(r.params[0][5:9]-1089.3)):
                ax1.fill_between(r.x,y,m*r.x/1000+b,color ='r')
                ax1top.fill_between(r.x,y-(m*r.x/1000+b),0,color ='r')
            else:
                ax1.plot(r.x,y,'k')
            
        tosave=append(tosave,[boundtodotpeak1],axis=0)
        tosave=append(tosave,[boundtodotpeak2],axis=0)
        tosave=transpose(tosave)
      #  savetxt('/home/chris/Ken/Cdenriched/'+z.name[-5]+'.csv',tosave,delimiter=',',header='x,baselined data,total fit'+(tosave.shape[1]-3)*',peak')



    
    xannotate = 1155
    yup = 100
    
    ax1.annotate('0eq', (xannotate,G[1160]+yup),horizontalalignment='right')
    ax1.annotate('25eq', (xannotate,F[1160]+yup),horizontalalignment='right')
    ax1.annotate('50eq', (xannotate,E[1160]+yup),horizontalalignment='right')
    ax1.annotate('80eq', (xannotate,D[1160]+yup),horizontalalignment='right')
    ax1.annotate('100eq', (xannotate,C[1160]+yup),horizontalalignment='right')
    ax1.annotate('200eq', (xannotate,B[1160]+yup),horizontalalignment='right')
    ax1.annotate('450eq', (xannotate,A[1160]+yup),horizontalalignment='right')
    ax1.set_xlim(1040,1160)
    ax1.set_ylim(-100,8000)
    
    
    
    
     
    
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
        tosave = array([r.x])
        tosave = append(tosave,[r.y_in],axis=0)        
        tosave = append(tosave,[r.y],axis=0)
        tosave=append(tosave,r.peaks,axis=0)
        tosave=transpose(tosave)
        savetxt('/home/chris/Ken/Stoichiometric/'+z.name[-5]+'.csv',tosave,delimiter=',',header='x,baselined data,total fit'+(tosave.shape[1]-3)*',peak')
        m = r.params[0][-2]
        b =r.params[0][-1]
        for p in range(len(r.peaks)):
            if p ==1:
                fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='g')
                if z is B:
                    ax2top.fill_between(r.x,r.peaks[p]-(m*r.x/1000+b),0,color ='g')
            elif p ==2:
                fill_between(r.x,r.peaks[p],m*r.x/1000+b,color ='b')
                if z is B:
                    ax2top.fill_between(r.x,r.peaks[p]-(m*r.x/1000+b),0,color ='b')
            
                    
            else:
                plot(r.x,r.peaks[p],'k')
        plot(r.x,r.y,'k')
    for z in [E] :
        r = fitspectrum(z, (1047,1100),'xGaussian', [500,1500,1500,
                                                         1064.76,  1078.21,  1085.099,
                                                        50,15,15,0,z[1100]])
        print r.params[0]
        tosave = array([r.x])
        tosave = append(tosave,[r.y_in],axis=0)        
        tosave = append(tosave,[r.y],axis=0)
        tosave=append(tosave,r.peaks,axis=0)
        tosave=transpose(tosave)
        savetxt('/home/chris/Ken/Stoichiometric/'+z.name[-5]+'.csv',tosave,delimiter=',',header='x,baselined data,total fit'+(tosave.shape[1]-3)*',peak')
        m = r.params[0][-2]
        b =r.params[0][-1]
        for p in range(len(r.peaks)):

            ax2.plot(r.x,r.peaks[p],'k')
        ax2.plot(r.x,r.y,'k')
    
        

    
    mbt.plot(color ='k',linestyle='--',ax=ax1top)
    cdmbt.plot(color = 'k',ax=ax1top)

    mbt.plot(color ='k',linestyle='--',ax=ax2top)
    cdmbt.plot(color = 'k',ax=ax2top)
  
    ax1top.set_ylim(0,1500)
    ax1top.set_xlim(1040,1160)
    ax2top.set_ylim(0,1700)
    ax2top.set_xlim(1040,1160)
    
    ax1top.set_yticks([])
    ax2top.set_yticks([])
    ax2.set_yticks([])
    ax1.set_yticks([])
    
    ax1top.set_xlabel('Raman Shift (cm$^{-1}$)')
    ax2top.set_xlabel('Raman Shift (cm$^{-1}$)')
    ax1.set_xlabel('Raman Shift (cm$^{-1}$)')
    ax2.set_xlabel('Raman Shift (cm$^{-1}$)')    
            
    return 0
    
    
    
def SH():
    
    
    
    os.chdir('/home/chris/Documents/DataWeiss/150728')
    cdmbt=copy.copy(CdMethylTPRef)
    mbt=copy.copy(MethylTPRef)
    cdmbt.autobaseline((1000,1200),order = 0)
    mbt.autobaseline((1000,1200),order = 0)

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
    for z in [A,B,C,D,E,F,G]:

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
        #z.to_csv('/home/chris/Dropbox/Ken/SHregion/'+z.name[-5]+'.csv')
        pass
 
    
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


    