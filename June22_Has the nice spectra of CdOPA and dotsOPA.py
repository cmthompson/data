# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:57:11 2015

@author: chris
"""
from numpy import *
from matplotlib.pyplot import *
from copy import deepcopy
from ramanTools.RamanSpectrum import *
from scipy.optimize import minimize
from matplotlib import gridspec

def June22():
    d = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_06.txt')   ###  pH 1
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_07.txt')   ### pH 5
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_08.txt')  ### pH 12
    c= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt')  ## dots
    d.autobaseline((283,1989),order = 3)
    d.autobaseline((284,408,577,701,826,1147,1380,1588,1849,1976),join='start',order = 9,specialoption='points')

    
    a.autobaseline((283,1989),order = 3)
    a.autobaseline((284,577,701,826,1147,1380,1588,1849,1976),join='start',order = 9,specialoption='points')
    
   
    b.autobaseline((283,1989),order = 3)
    b.autobaseline((284,577,701,826,1147,1380,1588,1849,1976),join='start',order = 9,specialoption='points')
    
   
    c= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt')  ## dots
    c.autobaseline((911,1196,1385,1515,1800),specialoption='points',order=7,join='start')
    c.autobaseline((280,600,690,826,861,911),specialoption='points', order = 5,join='end')
    
    c.autobaseline((281,630), order = 4, join = 'end')    
    c.autobaseline((1492,1800), order = 5, join = 'start')
    c[:]*=3
    
    c[:]+=5000
    b[:]+=4000
    a[:]+=1300
    c.plot()
    b.plot()
    a.plot()
    d.plot()
    ylim(0,10000)
    xlim(400,1800)
    
    legend(['dots','CdOPA pH12', 'pH 5', 'pH1'])

    tosave = transpose([d.index,d.values,a.values,b.values])
    savetxt('/home/chris/Desktop/emily/CdOPARaman.csv', tosave, header = 'pH1, pH5, pH12',delimiter = ',')
    c.to_csv('/home/chris/Desktop/emily/DotsOPARaman.csv')
    ylabel('Intensity (a.u.)')
    xlabel('Raman Shift (cm$^{-1}$)')
   # savefig('/home/chris/Dropbox/GroupmeetingJuly9_2015/dotsandrefs.png', dpi=256)
    return 0
    
def MBTSeries():
    clf()
    ax1 = gca()
    
    chdefarea=array([])
    thiolarea=array([])
    
    
    native= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt') ###### Native ligand only
    native[:]/=2
    #native=removespikes(native)
    native.autobaseline((911,1196,1385,1515,1800),join='start',specialoption='points',order=7)
    native.autobaseline((600,690,826,861,911),specialoption='points', order = 5,join='end')
   
    
    eightyfour= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150623/150623_7.txt') ###### 84 eq
    eightyfour[:]/=2
    #eightyfour=removespikes(eightyfour)
    eightyfour.autobaseline((764,838),order = 0,join='start')
    eightyfour.autobaseline((838,2000),order = 1,join='start')
    eightyfour.autobaseline((600,690,826,861,900,1196,1385,1515,1657),specialoption='points',order=7)
   # eightyfour.smooth()
    

    
    
    sixhundredforty= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150617/150617_01.txt')  ###### 640 eq MBT
    #sixhundredforty=removespikes(sixhundredforty)
    sixhundredforty.autobaseline((803,861),order=1, join='start')
    sixhundredforty.autobaseline((861,1254),order = 2, join='start')
    sixhundredforty.autobaseline((1254,1515),order = 4, join = 'start')
    sixhundredforty.autobaseline((1515,2000),order = 3, join = 'start')
    sixhundredforty.autobaseline((555,613,764,1141,1321,1565,1652),specialoption='points',order=6)
    #sixhundredforty.smooth()
    
    thirtytwo=RamanSpectrum('/home/chris/Dropbox/DataWeiss/150623/150623_10.txt')  ###### 32 equivalents
    #thirtytwo=removespikes(thirtytwo)
    thirtytwo.autobaseline((300,862),order = 1,join='start')
    thirtytwo.autobaseline((786,862),order = 0,join='start')
    thirtytwo.autobaseline((862,1425),order = 2,join='start')
    thirtytwo.autobaseline((1425,1439),order = 1,join='start')
    thirtytwo.autobaseline((1439,2000),order = 2,join='start')
    x = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150623/150623_12.txt')
    #x= removespikes(x)
    x.autobaseline((740,1441),order = 1)
    thirtytwo=add_RamanSpectra(thirtytwo,x)
    thirtytwo.autobaseline((740,764,1052,1141,1321,1425,1441),specialoption='points',order=2)
   # thirtytwo.smooth()
    
    
    oneseventynine=RamanSpectrum('/home/chris/Dropbox/DataWeiss/150623/150623_16.txt')  ###### 179 equivalents
   # oneseventynine=removespikes(oneseventynine)
    oneseventynine.autobaseline((300,793),order = 1,join='start')
    oneseventynine.autobaseline((793,862),order = 1,join='start')
    oneseventynine.autobaseline((862,1460),order = 1,join='start')
    oneseventynine.autobaseline((1460,1486),order = 1,join='start')
    oneseventynine.autobaseline((1486,2000),order = 1,join='start')
    oneseventynine.autobaseline((740,764,1052,1141,1321,1425,1441,1700),specialoption='points',order=2)
    #oneseventynine.smooth()
    
    mbt = CdMethylTPRef.copy()
    mbt[:]/=10
    
    lw = 2
    thirtytwo[:]+=200
    eightyfour[:]+=700
    oneseventynine[:]+=1000
    sixhundredforty[:]+=1550
    mbt[:]+=1950
    
    
    chtwistarea=array([native.calc_area((1285,1332)),thirtytwo.calc_area((1285,1332)),eightyfour.calc_area((1285,1332)),oneseventynine.calc_area((1285,1332)),sixhundredforty.calc_area((1285,1332)),mbt.calc_area((1285,1332))])
    chdefarea=array([native.calc_area((1413,1475)),thirtytwo.calc_area((1413,1475)),eightyfour.calc_area((1413,1475)),oneseventynine.calc_area((1413,1475)),sixhundredforty.calc_area((1413,1475)),mbt.calc_area((1413,1475))])
    thiolarea1=array([native.calc_area((1587,1611)),thirtytwo.calc_area((1587,1611)),eightyfour.calc_area((1587,1611)),oneseventynine.calc_area((1587,1611)),sixhundredforty.calc_area((1587,1611)),mbt.calc_area((1587,1611))])
    
    fits1 = list()
    fits2 = list()
    
    a = [thirtytwo, eightyfour, oneseventynine, sixhundredforty, mbt]
    a.reverse()
    for i in a:
        i.plot(linewidth=lw,axes=ax1)
    native.plot(linewidth = lw,axes=ax1)
    
    ax1.set_ylabel('Intensity (a.u.)')
    ax1.set_xlabel('Raman shift (cm$^{-1}$')
    legend(['solid', '640eq','179', '84eq','32eq','0'])
    ax1.set_xlim(500,1800)
    ax1.set_ylim(0,10000)
    

    
    ax2 = figure().add_subplot(111)
    for i in [thirtytwo, eightyfour, oneseventynine, sixhundredforty]:
        guess = [100,500,500,1065,1080,1085,7, 7,7,0,i[1100]]
        r = fitspectrum(i,(1050,1105), 'ThreeGaussian', guess)
        for p in r.peaks:
            ax1.plot(r.x, p,'k', linewidth = 2)
        fits1.append(r.areas[1]/r.areas[2]) 
    ax2.plot([32,84,179,640],fits1,'rs-', label='1')
    
    return None
    
    
def June29():
     a = loadtxt('/home/chris/Documents/DataWeiss/150629/after exchange.csv', delimiter = ',', unpack = True,skiprows=1)
     a=a[:,:450]
     for i in a[2:5]:
         i-=min(i)
         i/=max(i)
         
     plot(a[0],a[1],label = 'glass slide')
     plot(a[0],a[2],label = 'dots in water (control)')
     plot(a[0],a[3],label = 'dots NiCl2 anhydrous')
     plot(a[0],a[4],label = 'dots NiCl2 hexahydrate')
     legend()
     
#     figure()
#     a = loadtxt('/home/chris/Documents/DataWeiss/150629/Solutions of NiCl2 after exchange.csv', delimiter = ',', unpack = True,skiprows=1)
#    # a=a[:,:450]
#     for i in a[1:3]:
#         i-=min(i)
#         i/=max(i)
#         
#     plot(a[0],a[1],label = 'NiCl anhydrous')
#     plot(a[0],a[2],label = 'NiCl2 hexahydrate')
#
#     legend()
     return 0
    
    
def June30():
    
    d = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_06.txt')   ###  pH 1
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_07.txt')   ### pH 5
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_08.txt')  ### pH 12
    c= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt')  ## dots
    d.autobaseline((283,1989),order = 3)
    a.autobaseline((283,1989),order = 3)
    b.autobaseline((283,1989),order = 3)
    
    c.autobaseline((600,690,826,861,900,1196,1385,1515,1657),specialoption='points',order=7)
    
   #c[:]+=4500
    #b[:]+=3000
    #a[:]+=1500
    c.plot()
    b.plot()
    a.plot()
    #d.plot()
    
    legend(['dots','CdOPA pH12', 'pH 5', 'pH1'])
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150630/150623_2.txt')
    
    a.autobaseline((700,1500), order  = 0)
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150630/150623_3.txt')
    b.autobaseline((300,791), order =2)
    b.autobaseline((791,858),order = 0,join='start')
    b.autobaseline((858,2000),order =1,join='start')
    b.autobaseline((400,700,954,1495,1700),specialoption='points',order=3)
    #a.plot()
    c = add_RamanSpectra(a,b)
    
    b.plot()
    OPARef.plot()
    return 0   
    
def June30MetalExchange():
    a = loadtxt('/home/chris/Documents/DataWeiss/150630/Metal Exchange exp after one day 150630.csv', delimiter = ',', unpack = True,skiprows = 2)
    a = a[:,:850]
    a[0]= 10**7/a[0]
    for i in range(1,a.shape[0]):
        plot(a[0],a[i])
    legend(['Control - No liquid',	'Acidonly ph5',	'Cu exchanged',	'Mn exchanged',	'Nickel exchanged',	'lead exchanged'])
    return 0
    
    
#def June30Raman():
#    MBTSeries()
#    a = sixhundredforty= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150630/SampleA640eq.txt')  ###### 640 eq MBT retaken on June30
#    sixhundredforty=removespikes(sixhundredforty)
#    sixhundredforty.autobaseline((109,500),order=1, join='start')
#    sixhundredforty.autobaseline((500,723),order=1, join='start')
#    sixhundredforty.autobaseline((723,795),order=1, join='start')
#    sixhundredforty.autobaseline((795,1930),order = 2, join='start')
#    sixhundredforty.autobaseline((200,555,613,764,1052,1141,1321,1565,1700),specialoption='points',order=6)
#    sixhundredforty.smooth()
#    
#    sixhundredforty[:]+=1100
#    
#    
#    sixhundredforty.plot(linewidth=2)
#    
#    a = thirtytwo= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150630/150630_sampleD_128x10s.txt')  ###### 640 eq MBT retaken on June30
#    thirtytwo=removespikes(thirtytwo)
#    thirtytwo.autobaseline((109,500),order=1, join='start')
#    thirtytwo.autobaseline((500,723),order=1, join='start')
#    thirtytwo.autobaseline((723,795),order=1, join='start')
#    thirtytwo.autobaseline((795,1930),order = 2, join='start')
#    thirtytwo.autobaseline((200,555,613,764,1052,1141,1321,1565,1700),specialoption='points',order=6)
#    thirtytwo.smooth()
#    
#    thirtytwo[:]+=100
#    
#    
#    thirtytwo.plot(linewidth=2)
#    
#    return 0
    
def July1():
    
   
    
    
    ratiolist = list()
    native= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt') ###### Native ligand only
    native[:]/=2
    native=removespikes(native)
    native.autobaseline((600,690,826,861,900,1196,1385,1515,1657),specialoption='points',order=7)
    native.smooth()
    
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files1.txt') 
    b= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files2.txt')
    c= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files3.txt')
    d= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files4.txt')
    e= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files5.txt')
    
    #c = removespikes(c)
    #c = removespikes(c)
    
    
    correct = zeros(e.values.shape)
    for z in [a,b,c,d,e]:
        y = deepcopy(z)
        y.smooth()
        y.smooth()
        y.smooth()
        correct+= y/z
    correct/=5 

    ax1=figure().add_subplot(111)
    mbt = CdMethylTPRef.copy()
    mbt[:]/=10
    

    for z in [a,b,c,d,e]:

        z[:]*=correct
        z=removespikes(z)
        z.autobaseline((109,500),order=3, join='start')
        z.autobaseline((500,725),order=2, join='start')
        z.autobaseline((725,795),order=1, join='start')
        z.autobaseline((795,1363),order=2, join='start')
        z.autobaseline((1363,1430),order = 1, join = 'start')
        z.autobaseline((1430,1930),order = 4, join='start')
        z.autobaseline((200,555,613,764,1141,1321,1565,1700,1920),specialoption='points',order=7)
        z.smooth()

    
    mbt[:]+=3000
    a[:]+=1000
    b[:]+=600
    c[:]+=400
    d[:]+=200
    e[:]-=200
    
    native-=500
    
    lw = 2
    mbt.plot(linewidth = lw)
    guess = [50,100,100,1065,1080,1085,7,7,7,0,z[1110]] 
    for z in [a,b,c,d,e]:
        print z.name

        
        r = fitspectrum(z,(1050,1110),'ThreeGaussian',[50,100,500,1065,1078,1085,15,15,15,0,z[1110]] ) #r = fitspectrum(z,(1070,1110),'TwoGaussian',guess )
        
        ratio = r.areas[1]/r.areas[2]
        
        if z is a:
            r = fitspectrum(z,(1070,1110),'TwoGaussian',[100,1000,1078,1085,15,15,0,z[1110]] ) #r = fitspectrum(z,(1070,1110),'TwoGaussian',guess )
            ratio = r.areas[0]/r.areas[1]
        ratiolist.append(ratio)
        z.plot(linewidth = lw)
        print r.params[0]
        for p in r.peaks:
            ax1.plot(r.x,p,color = 'k',linewidth = 2)
        plot(r.x,r.y, color = 'k', linewidth = lw)
          
   


    native.plot(linewidth = lw)
    
    xlim(555,1700)
    ylim(-500,3000)
    
    #legend(['mbt solid','2035eq','713eq','502eq','80eq','58eq','native'])
    
    ax2=figure().add_subplot(111)
    print ratiolist
    ax2.plot([2035,713,502,80,58],ratiolist,'rs-')

    
    
    return 0
    
def July1phosphonicacidtreated():
    
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_07.txt')   ### pH 5
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150622/150622_08.txt')  ### pH 12
    c= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt')  ## dots
   
    a.autobaseline((283,1989),order = 3)
    b.autobaseline((283,1989),order = 3)
    
    c.autobaseline((600,690,826,861,900,1196,1385,1515,1657),specialoption='points',order=7)
    
    c[:]+=4500
    b[:]+=3000
    a[:]+=1500
    c.plot()
    b.plot()
    a.plot()

    c= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150701/files6phosphonate.txt'  )
    c.autobaseline((283,1989),order = 3)
    c[:]*=3
    c.plot()
    
    
    a = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150630/150623_2.txt')
    
    a.autobaseline((700,1500), order  = 0)
    b = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150630/150623_3.txt')
    b.autobaseline((300,791), order =2)
    b.autobaseline((791,858),order = 0,join='start')
    b.autobaseline((858,2000),order =1,join='start')
    b.autobaseline((400,700,954,1495,1700),specialoption='points',order=3)
    
    c = add_RamanSpectra(a,b)
    
    x = RamanSpectrum('/home/chris/Dropbox/DataWeiss/150623/150623_4.txt')
    b.plot()
    OPARef.plot()
    x.plot()
    legend(['dots','CdOPA pH12', 'pH 5','phosJuuly1', 'phosJune23', 'OPA ref','june23'])
    return 0   
    
def OPAMBTExchange():
    figure()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150707/150707_02.txt')
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150707/150707_03.txt')
    a[:]*=5
    b[:]*=5
    subplot(121)
    a.plot()
    CdMethylTPRef.plot()
    legend(['new','old'])
    r = fitspectrum(a,(1070,1110),'OneGaussian',[25000,1088,10,0,0])
    plot(r.x,r.y, 'k',linewidth = 2)
    
    r = fitspectrum(CdMethylTPRef,(1070,1110),'OneGaussian',[25000,1088,10,0,0])
    plot(r.x,r.y,'r' ,linewidth = 2)    
    
    subplot(122)
    b.plot()
    CdMeOTPRef.plot()
    legend(['new','old'])
    
    r = fitspectrum(b,(1070,1110),'OneGaussian',[60000,1088,10,0,0])
    plot(r.x,r.y,'k',linewidth = 2)
    
    r = fitspectrum(CdMeOTPRef,(1070,1110),'OneGaussian',[6000,1088,10,0,0])
    plot(r.x,r.y,'r', linewidth = 2)    
    
    figure()
    June22()
    c= RamanSpectrum('/home/chris/Documents/DataWeiss/150707/150707_05.txt')
    
    c.autobaseline((200,2000),order= 4)
    c[:]+=4000
    c.plot()
    
    
    figure()
    d= RamanSpectrum('/home/chris/Documents/DataWeiss/150707/150707_06.txt')
    d.autobaseline((200,2000),order= 4)
    d[:]*=10
    d.plot()
    a.plot()
    legend(['exchanged', 'reference'])
    return 0
    
def CdMBTinDMF():
    clf()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150709/150709_01.txt')#### DMFonly
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150709/150709_02.txt')#### 510mg DMF with 200 mgCdMBT
    a.autobaseline((523,935,1336,1780),order = 3,specialoption='points')
    b.autobaseline((523,935,1336,1780),order = 3,specialoption='points')
    a[:]*=4720
    a[:]/=6256
   
    c = RamanSpectrum(b-a)
    c.plot()
    r = fitspectrum(c,(1070,1105),'xVoigt',[10000,1088,15,6,0,0])
    plot(r.x,r.y,'s-',linewidth=2)
    for i in r.peaks:
        plot(r.x,i)
    print r.areas
    print r.params[0][2:4]
    CdMethylTPRef.plot()
        
    
#    def difference(c): return sum((b[200:1700]-c*a[200:1700])**2)
#    r = minimize(difference,[1])
    
    return r.params
    
def DMFWash():
    clf()
    ax1 = gca()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150709/150709_03.txt')#### sample A washed with DMF
    
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150709/150709_04.txt')# sa,[;e B washed with DMF]
    
    c = RamanSpectrum('/home/chris/Documents/DataWeiss/150709/150709_05.txt')##sample C washed with DMF
    d = RamanSpectrum('/home/chris/Documents/DataWeiss/150709/150709_06.txt')##sample C washed with DMF using 50x close up objective
    e= RamanSpectrum('/home/chris/Documents/DataWeiss/150709/150709_07.txt')##sample E washed with DMF
 
    for z in [a,b,c]:
        z.autobaseline((200,725),order=3, join='start')    
        z.autobaseline((725,800),order=0, join='start')
        z.autobaseline((800,1427),order=2, join='start')
        z.autobaseline((1427,1435),order=0, join='start')
        z.autobaseline((1435,2000),order=0, join='start')   
        z[:]-=z[1700]
        z.plot()
        
    d.autobaseline((520,1250),order = 3)
    e.autobaseline((520,1250),order = 3)
      
    d.plot()
    e.plot()
    xlim(900,1200)
    
    ax1=figure().add_subplot(111)

    
    ratiolist = list()
#    native= RamanSpectrum('/home/chris/Dropbox/DataWeiss/150612/150612_01_CdSe.txt') ###### Native ligand only
#    native[:]/=2
#    native=removespikes(native)
#    native.autobaseline((600,690,826,861,900,1196,1385,1515,1657),specialoption='points',order=7)
#    native.smooth()
#    
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files1.txt') 
    b= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files2.txt')
    c_unwashed= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files3.txt')
    d_unwashed= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files4.txt')
    e_unwashed= RamanSpectrum('/home/chris/Documents/DataWeiss/150701/files5.txt')
    

    
    
    correct = zeros(e_unwashed.values.shape)
    for z in [a,b,c_unwashed,d_unwashed,e_unwashed]:
        y = deepcopy(z)
        y.smooth()
        y.smooth()
        y.smooth()
        correct+= y/z
    correct/=5 

    
    

    for z in [c_unwashed,e_unwashed]:

        z[:]*=correct
        z=removespikes(z)
        z.autobaseline((109,500),order=3, join='start')
        z.autobaseline((500,725),order=2, join='start')
        z.autobaseline((725,795),order=1, join='start')
        z.autobaseline((795,1363),order=2, join='start')
        z.autobaseline((1363,1430),order = 1, join = 'start')
        z.autobaseline((1430,1930),order = 4, join='start')
        z.autobaseline((200,555,613,764,1141,1321,1565,1700,1920),specialoption='points',order=7)
        z.smooth()

    lw = 2
    c_unwashed[:]*=3

    guess = [300,1500,1078,1085,15,15,0,0]

    r = fitspectrum(d,(1070,1110),'TwoGaussian',guess ) #r = fitspectrum(z,(1070,1110),'TwoGaussian',guess )
    ratio = r.areas[0]/r.areas[1]
    ratiolist.append(ratio)
    d.plot(linewidth = lw)
    print r.params[0]
    for p in r.peaks:
        ax1.plot(r.x,p,color = 'k',linewidth = 2)
    plot(r.x,r.y, color = 'k', linewidth = lw)
    
    r = fitspectrum(c_unwashed,(1070,1110),'TwoGaussian', guess) #r = fitspectrum(z,(1070,1110),'TwoGaussian',guess )
    ratio = r.areas[0]/r.areas[1]
    ratiolist.append(ratio)
    c_unwashed.plot(linewidth = lw)
    print r.params[0]
    for p in r.peaks:
        ax1.plot(r.x,p,color = 'k',linewidth = 2)
    plot(r.x,r.y, color = 'k', linewidth = lw)
    
    
    e[:]+=2000
    e_unwashed[:]+=2000
    e.smooth()
    guess = [300,1500,1078,1085,15,15,0,2000]

    r = fitspectrum(e,(1050,1110),'ThreeGaussian',  [100,300,1500,1065,1078,1085,15,15,15,0,2000] ) #r = fitspectrum(z,(1070,1110),'TwoGaussian',guess )
    ratio = r.areas[1]/r.areas[2]
    ratiolist.append(ratio)
    e.plot(linewidth = lw)
    print r.params[0]
    for p in r.peaks:
        ax1.plot(r.x,p,color = 'k',linewidth = 2)
    plot(r.x,r.y, color = 'k', linewidth = lw)
    
    r = fitspectrum(e_unwashed,(1070,1110),'TwoGaussian',guess) #r = fitspectrum(z,(1070,1110),'TwoGaussian',guess )
    ratio = r.areas[0]/r.areas[1]
    ratiolist.append(ratio)
    e_unwashed.plot(linewidth = lw)
    print r.params[0]
    for p in r.peaks:
        ax1.plot(r.x,p,color = 'k',linewidth = 2)
    plot(r.x,r.y, color = 'k', linewidth = lw)
        
    
    legend(['c-washed1', 'c-unwashed','e-washed', 'eunwashed'])
    
    print ratiolist
    return 0
    
    
def July10():
    A = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_01.txt')
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_02.txt')
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_03.txt')
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_04.txt')
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_05.txt')
    Awash = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_06.txt')
    Bwash = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_07.txt')
    Cwash = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_08.txt')
    Dwash = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_09.txt')
    
  
    for z in [A,B,C,D,E,Awash,Bwash,Cwash,Dwash]:
        
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
        z[:]-=z[1667]
        
    
        
    A[:]+=20000
    B[:]+=15000
    C[:]+=10000
    D[:]+=5000
    E[:]+=0
    
    Awash[:]+=20000
    Bwash[:]+=15000
    Cwash[:]+=10000
    Dwash[:]+=5000
    
    
    
    A.plot(color = 'r')
    B.plot(color = 'r')
    C.plot(color = 'r')
    D.plot(color = 'r')
    E.plot(color = 'r')
    
    
    Awash.plot(color = 'k')
    Bwash.plot(color = 'k')
    Cwash.plot(color = 'k')
    Dwash.plot(color = 'k')
    
    ratiolist = list()
    for z in [A,B,C,D,E]:
        pass
        #r = fitspectrum(z,(1000,1130),'xGaussian', [100,100,100,100])
    
    ylim(0,30000)
   
    return 0
    
def July10plusRef():
    fig = figure(figsize=(12, 6)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
    ax1 = subplot(gs[0])
   
    A = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_01.txt')
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_02.txt')
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_03.txt')
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_04.txt')
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150710/150710_05.txt')
   
    
  
    for z in [A,B,C,D,E]:
        
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
        z[:]-=z[1667]
        
    A[:]+=7000
    B[:]+=4500
    C[:]+=3500
    D[:]+=2000
    E[:]+=1100
    
    (ODPARef/6+8000).plot()
    (MethylTPRef/20+8000).plot()
    (CdMethylTPRef/2+6000).plot()
    A.plot()
    B.plot()
    C.plot()
    D.plot()
    E.plot()

    ax1.set_ylim(0,16000)
    ax1.set_xlim(500,1800)
    
    ax2 = subplot(gs[1])
    
    A.plot()
    B.plot()
    C.plot()
    D.plot()
    E.plot()
    ax2.set_xlim(1050,1110)
    ax2.set_ylim(0,10000)
    ax2.annotate('*', (1079,9500))
    ax2.annotate('*', (1086,9500))
    ax2.annotate('0eq', (1100,620))
    ax2.annotate('50eq', (1100,1700))
    ax2.annotate('100eq', (1100,3300))
    ax2.annotate('500eq', (1100,4600))
    ax2.annotate('1000eq', (1100,7500))
    
    for z in [A,B,C,D,E] :
        r = fitspectrum(z, (1047,1100),'xGaussian', [500,1500,1500,
                                                         1064.76,  1078.21,  1085.099,
                                                        50,15,15,0,z[1100]])
        print r.params[0][6:9]
        print 'ratio of A',  r.params[0][1]/r.params[0][2]
        
        for p in r.peaks:
            plot(r.x,p,'k')
        plot(r.x,r.y,'k')
    
        
   
    return 0
def July10UVVis():
    a = loadtxt('/home/chris/Documents/DataWeiss/150710/stoichdots after MBT treatement.csv',delimiter = ',', skiprows = 2, unpack = True)
    x1 = argmin(abs(a[0]-505))
    x2 = argmin(abs(a[0]-515))
    xfoot = argmin(abs(a[0]-560))
    for i in a[4:8]:
        i = SGsmooth(a[0],i)

        i-=i[xfoot]
        ymax= max(i[x2:x1])
        
        i/=ymax
        plot(a[0],i)
        
    i = a[3]
    i-=i[xfoot]
    ymax= max(i[x2:x1])
    
    i/=ymax
    plot(a[0],i)
    legend(['1000','500','100','50','0'])
    return 0
    
def July11():
    fig = figure(figsize=(12, 6)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
    ax1 = subplot(gs[0])
    
    A = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_01.txt')
    B = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_02.txt')
    C = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_03.txt')
    D = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_04.txt')
    E = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_05.txt')
    Awash = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_06.txt')
    Bwash = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_07.txt')
    Cwash = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_08.txt')
    Dwash = RamanSpectrum('/home/chris/Documents/DataWeiss/150711/150711_09.txt')
    
  
    for z in [A,B,C,D,E,Awash,Bwash,Cwash,Dwash]:
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
        
    A[:]+=15000
    B[:]+=6000
    C[:]+=4000
    D[:]+=2000
    E[:]+=0
    
    Awash[:]+=8000
    Bwash[:]+=6000
    Cwash[:]+=4000
    Dwash[:]+=2000
       
    

    
    
    
    A.plot(color = 'r')
    B.plot(color = 'r')
    C.plot(color = 'r')
    D.plot(color = 'r')
    E.plot(color = 'r')
    
    
#    Awash.plot(color = 'k')
#    Bwash.plot(color = 'k')
#    Cwash.plot(color = 'k')
#    Dwash.plot(color = 'k')
    
    ylim(0,30000)
    
    ax2 = subplot(gs[1])
    
    A.plot(color = 'k')
    B.plot(color = 'k')
    C.plot(color = 'k')
    D.plot(color = 'k')
    E.plot(color = 'k')
    
    ax2.set_xlim(1050,1110)
    ax2.set_ylim(0,30000)
    #ax2.annotate('*', (1079,9500))
    #ax2.annotate('*', (1086,9500))
    ax2.annotate('0eq', (1100,240))
    ax2.annotate('50eq', (1100,2200))
    ax2.annotate('100eq', (1100,4400))
    ax2.annotate('500eq', (1100,6400))
    ax2.annotate('1000eq', (1100,12400))
    
    

    
    
    rlist = list()
#    for z in [C,D]:
#        r = fitspectrum(z, (1000,1160),'xGaussian', [200,200,500,500,500,100,100,500,
#                                                        1032.25,  1062.17,1070,1078,1086,1116.75, 1129.13, 1133.6,
#                                                        50,15,50,15,15,15,15,15,0,z[1160]])
##        plot(r.x,r.y,color = 'b', linewidth = 2)
#            for i in r.peaks:
#            plot(r.x,i)
    g0a = 13.8   #### from july10data average of four spectra widhts
    g0b = 21.3  ### from July10data       
    Aratio = 0.875

        
    
    def func(x,A0, A1,A2,A3,A4,A5,A6,A7,w1,w2,w3,w4,w5,w6,w7,G1,G2,G3,G4,G5,G6,G7,m,b): 
        return m*x/1000+b + A0*Aratio*exp(-(x-1078.1)**2/g0a)+A0*exp(-(x-1085.5)**2/g0b)+A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2)  +A3*exp(-(x-w3)**2/G3)  + A4*exp(-(x-w4)**2/G4)  + A5*exp(-(x-w5)**2/G5)  + A6*exp(-(x-w6)**2/G6)+ A7*exp(-(x-w7)**2/G7)
    
    dotbound_polymerbound_ratio=list()    
    
    
    for z in [A,B]:
        r = fitspectrum(z, (1000,1160),'Custom', [300, 
                                                        387,      740,       568,     7000,  241  ,   334,     1309,
                                                        1032.25,  1061,  1062.17,1089, 1116.75, 1129.13, 1133.6,
                                                        150,25,1000,15,15,15,15,0,z[1160]],function=func)
        print r.params[0]
        A0 = r.params[0][0]
        m = r.params[0][-2]
        b =r.params[0][-1]
        boundtodotpeaks = A0*Aratio*exp(-(r.x-1078.1)**2/g0a)+A0*exp(-(r.x-1085.5)**2/g0b)+m*r.x/1000+b
        
        
        plot(r.x,r.y,color = 'b', linewidth = 2)
        
        peaks = list()
        areas = list()
        
        for i in range(7):
            A = r.params[0][1:][i]
            x0 = r.params[0][1:][i+7]
            G = r.params[0][1:][i+2*7]
            x= r.x
                
            
            y= A*exp(-(x-x0)**2/G)+r.params[0][1:][-2]/1000*x+r.params[0][1:][-1]
            ar1 = A*numpy.sqrt(pi*G)
            areas.append(ar1)
            plot(r.x,y,'k')
            if i ==3:
                fill_between(r.x,y,m*r.x/1000+b,color = 'b')
        ratioofdotboundareastopolymerboundarea = r.params[0][0]*numpy.sqrt(pi)*(Aratio*numpy.sqrt(13.8)+numpy.sqrt(21.3))/areas[3]
        dotbound_polymerbound_ratio.append(ratioofdotboundareastopolymerboundarea)
        
        fill_between(r.x,boundtodotpeaks,m*r.x/1000+b,color ='r')
    def func(x,A0, A1,A2,A3,A4,A5,A6,w1,w2,w3,w4,w5,w6,G1,G2,G3,G4,G5,G6,m,b): 
        return m*x/1000+b + A0*Aratio*exp(-(x-1078.1)**2/g0a)+A0*exp(-(x-1085.5)**2/g0b)+A1*exp(-(x-w1)**2/G1) +A2*exp(-(x-w2)**2/G2)  +A3*exp(-(x-w3)**2/G3)  + A4*exp(-(x-w4)**2/G4)  + A5*exp(-(x-w5)**2/G5)  + A6*exp(-(x-w6)**2/G6)
       
    for z in [C,D]:
        r = fitspectrum(z, (1000,1160),'Custom', [600, 
                                                        300,      740,       568,      241  ,   334,     1309,
                                                        1032.25,  1061,  1062.17, 1116.75, 1129.13, 1133.6,
                                                        150,25,600,15,15,15,0,z[1160]],function=func)
        print r.params[0]
        A0 = r.params[0][0]
        m = r.params[0][-2]
        b =r.params[0][-1]
        boundtodotpeaks = A0*Aratio*exp(-(r.x-1078.1)**2/g0a)+A0*exp(-(r.x-1085.5)**2/g0b)+m*r.x/1000+b
        plot(r.x,boundtodotpeaks)
        
        plot(r.x,r.y,color = 'k', linewidth = 2)
        
        peaks = list()
        areas = list()
        numpeaks = 6
        for i in range(numpeaks):
            A = r.params[0][1:][i]
            x0 = r.params[0][1:][i+numpeaks]
            G = r.params[0][1:][i+2*numpeaks]
            x= r.x
                
            
            y= A*exp(-(x-x0)**2/G)+r.params[0][1:][-2]/1000*x+r.params[0][1:][-1]
            ar1 = A*numpy.sqrt(pi*G)
            areas.append(ar1)
            plot(r.x,y,'k')
        
        fill_between(r.x,boundtodotpeaks,m*r.x/1000+b,color = 'r')
        ratioofdotboundareastopolymerboundarea = 100#r.params[0][0]*numpy.sqrt(pi)*(Aratio*numpy.sqrt(13.8)+numpy.sqrt(21.3))/areas[3]
        dotbound_polymerbound_ratio.append(ratioofdotboundareastopolymerboundarea)
    print dotbound_polymerbound_ratio
    figure()
    plot([1000,500,100,50],dotbound_polymerbound_ratio)
    return 0 
        
    for z in [E]:
        r = fitspectrum(z, (1000,1160),'xGaussian', [200,200,500,500,500,2000,100,100,500,500,
                                                        1032.25,  1050,   1062.17,1064.76,  1078.21,  1085.099,1090,1116.75, 1129.13, 1145,
                                                        50,50,15,15,15,15,15,15,15,15,0,z[1160]])
#
#        for p in r.peaks:
#            plot(r.x,p,color='b')
        plot(r.x,r.y,color = 'b', linewidth = 2)
    ax2.set_xlim(1000,1200)
    ax2.set_ylim(0,4000)       
#    figure()        
#    plot([50,100,500,1000],rlist)
    
   
    return 0    
    

def July13UVVis():
    subplot(121)
    annotate('Cd-enriched', (0.2,0.9), xycoords= 'axes fraction',fontsize = 24)
    a = loadtxt('/home/chris/Documents/DataWeiss/150713/150713-cdse enriched dots thiophenol treated.csv',delimiter = ',', skiprows =1, unpack = True)
    x1 = argmin(abs(a[0]-535))
    x2 = argmin(abs(a[0]-545))
    xfoot = argmin(abs(a[0]-620))
    for i in a[1:]:
        i = SGsmooth(a[0],i)

        i-=i[xfoot]
        ymax= max(i[x2:x1])
        
        i/=ymax
        plot(a[0],i)
        
    i = a[3]
    i-=i[xfoot]
    ymax= max(i[x2:x1])
    
    i/=ymax
    plot(a[0],i)
    legend(['500eq','250eq','100eq','50eq','25eq','0eq'])
    subplot(122)
    annotate('Stoichiometric', (0.2,0.9), xycoords = 'axes fraction',fontsize = 24)
    a = loadtxt('/home/chris/Documents/DataWeiss/150713/150713-cdse stoichiometric dots thiophenol treated.csv',delimiter = ',', skiprows =1, unpack = True)
    x1 = argmin(abs(a[0]-505))
    x2 = argmin(abs(a[0]-515))
    xfoot = argmin(abs(a[0]-560))
    for i in a[1:]:
        i = SGsmooth(a[0],i)

        i-=i[xfoot]
        ymax= max(i[x2:x1])
        
        i/=ymax
        plot(a[0],i)
        
    i = a[3]
    i-=i[xfoot]
    ymax= max(i[x2:x1])
    
    i/=ymax
    plot(a[0],i)
    legend(['500eq','250eq','100eq','50eq','25eq','25eq retake'])
    return 0