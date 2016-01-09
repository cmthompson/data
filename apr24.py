# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 19:39:13 2015

@author: chris
"""
from pandas import concat
import copy
from ramanTools.RamanSpectrum import *
def Apr24():
    os.chdir('/home/chris/Documents/DataWeiss/150424')
    a  = RamanSpectrum('150424_07.txt')
    b = RamanSpectrum('150424_08.txt')
    c  = RamanSpectrum('150424_09.txt')
    d  = RamanSpectrum('150424_10.txt')
    ax1 = gca()
    for i in (a,b,c,d):
        
        i = removespikes(i)
        
        i=smooth(i, window_len=9, window = 'SG')
        
        i = autobaseline(i,(112, 321),order = 2)
       
        i = autobaseline(i,(321, 766),order = 0,join='start')
        
        i = autobaseline(i,(766, 838),order = 0,join='start')
        i = autobaseline(i,(838, 1405),order = 2,join='start')
        i = autobaseline(i,(1405,1466),order = 0,join='start')
        
        i = autobaseline(i,(1466, 1974),order = 2,join='start')
        
        print calc_noise(i,(1700,1800))
        i.plot(marker = 'o',markersize = 1)
    for i in range(len(ax1.lines)):
        ax1.lines[i].set_ydata(ax1.lines[i].get_ydata()+i*200)
    
   
   
    (MeOTPRef/400+600).plot()
    
    legend(['0 eq', '50 eq', '100 eq', '200 eq'])

def Apr24b():
    
    os.chdir('/home/chris/Documents/DataWeiss/150424')
    a  = RamanSpectrum('150424_01.txt')
    b = RamanSpectrum('150424_04.txt')
    c  = RamanSpectrum('150424_10.txt')
    
    ax1 = gca()
    for i in (a,b,c):
        
        i = removespikes(i)
       
        #i = SPIDcorrect785(i)
        #i=smooth(i, window_len=9, window = 'SG')
        i = autobaseline(i,(112, 322),order = 2)
        i = autobaseline(i,(322, 767),order = 0,join='start')
        i = autobaseline(i,(767, 838),order = 0,join='start')
        i = autobaseline(i,(838, 1405),order = 2,join='start')
        i = autobaseline(i,(1405,1466),order = 0,join='start')
        
        i = autobaseline(i,(1466, 1974),order = 2,join='start')
        
        
        i.plot(marker = 'o',markersize = 1)
    for i in range(len(ax1.lines)):
        ax1.lines[i].set_ydata(ax1.lines[i].get_ydata()+i*200)
   
   
    legend(['chlorothiophenol', 'methylbenzenethiol', 'methoxythiophenol'])
    return 0
    
def Apr24UVVis():
    a = loadtxt('/home/chris/Documents/DataWeiss/150424/titration and Chloro and methyl.csv',delimiter = ',', unpack=True,skiprows=2)
    
    lstyle= '-'
    fig = gcf()
    ax1=fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    print a[0][230]
    for i in [1,2,3,4,6,8]:
        a[i] = SGsmooth(a[0],a[i])
        peakmax= argmax(a[i,200:240])+200
        print a[0,peakmax]
        a[i]-=min(a[i])
        a[i]/=a[i,peakmax]
        
        ax1.plot(a[0],a[i],lstyle, linewidth=3)
        ax2.plot(a[0],a[i],lstyle, linewidth=3)
    
    
    ax1.legend(['0 eq MTP', '50 eq MTP','100 eq MTP', '200 eq MTP', 'Cl',  'MBT'])
        
    ax1.set_ylim(0,2)
    ax1.set_xlim(450,800)
    ax1.set_title('570 nm dots with MTP, MBT, and ClTP')
    ax1.set_xlabel('wavelength (nm)')
    ax1.set_ylabel('absorbance (a.u.)')
    
    
    ax2.set_ylim(0.97,1.02)
    ax2.set_xlim(562,583)
    ax2.set_title('570 nm dots with MTP, MBT, and ClTP')
    ax2.set_xlabel('wavelength (nm)')
    ax2.set_ylabel('absorbance (a.u.)')

    return 0
    
def Apr30():
        
    os.chdir('/home/chris/Documents/DataWeiss/150424')
    a  = RamanSpectrum('150424_01.txt')
    b = RamanSpectrum('150424_04.txt')
    c  = RamanSpectrum('150424_10.txt')
    
    ax1 = gca()
    for i in (a,b,c):
        
        i = removespikes(i)
       
        #i = SPIDcorrect785(i)
        #i=smooth(i, window_len=9, window = 'SG')
        i.autobaseline((68, 322),order = 4)
        i.autobaseline((322, 767),order = 0,join='start')
        i.autobaseline((767, 838),order = 0,join='start')
        i.autobaseline((838, 1405),order = 2,join='start')
        i.autobaseline((1405,1466),order = 0,join='start')
        i.autobaseline((1466, 1974),order = 2,join='start')
        
        
        i.plot(marker = 'o',markersize = 1)
    
        
    os.chdir('/home/chris/Documents/DataWeiss/150430')
    a  = RamanSpectrum('150430_01.txt')
    b = RamanSpectrum('150430_03.txt')

    
    ax1 = gca()
    for i in (a,b):
        
        i = removespikes(i)

        i.autobaseline((68, 176),order = 2)
        i = autobaseline(i,(176, 767),order = 0,join='start')
        i = autobaseline(i,(767, 838),order = 0,join='start')
        i = autobaseline(i,(838, 1405),order = 2,join='start')
        i = autobaseline(i,(1405,1466),order = 0,join='start')
        
        i = autobaseline(i,(1466, 1980),order = 2,join='start')
        
        
        i.plot(marker = 'o',markersize = 1)
    
    ax1.lines[0].set_ydata(ax1.lines[0].get_ydata()+200)
    ax1.lines[1].set_ydata(ax1.lines[1].get_ydata()+500)
    ax1.lines[2].set_ydata(ax1.lines[2].get_ydata()/2+800)
    ax1.lines[3].set_ydata(ax1.lines[3].get_ydata()+1100)
    ax1.lines[4].set_ydata(ax1.lines[4].get_ydata()+1600)
    #ax1.lines[5].set_ydata(ax1.lines[5].get_ydata()+1000)
    ylim(0,1800)
    xlim(68,1980)
    vlines(1100,0,7500)
   
    xlabel('Raman Shift cm$^{-1}$')
    ylabel('Intensity')
    legend(['chlorothiophenol', 'methylbenzenethiol', 'methoxythiophenol','bromothiophenol', 'fluorothiophenol'])
    
    
    figure()
    ClTPRef.plot()
    (MethylTPRef/100).plot()
    (MeOTPRef/100).plot()
    BrTPRef.plot()
    FTPRef.plot()
    ax1=gca()
    ax1.lines[0].set_ydata(ax1.lines[0].get_ydata())
    ax1.lines[1].set_ydata(ax1.lines[1].get_ydata()+1500)
    ax1.lines[2].set_ydata(ax1.lines[2].get_ydata()+3000)
    ax1.lines[3].set_ydata(ax1.lines[3].get_ydata()+4500)
    ax1.lines[4].set_ydata(ax1.lines[4].get_ydata()+6000)
    ylim(0,7500)
    xlim(68,1980)
    vlines(1100,0,7500)
    xlabel('Raman Shift cm$^{-1}$')
    ylabel('Intensity')
    legend(['chlorothiophenol', 'methylbenzenethiol', 'methoxythiophenol','bromothiophenol', 'fluorothiophenol'])
    
    return 0
   


def thiophenolfits():
    
        
    os.chdir('/home/chris/Documents/DataWeiss/150424')
    a  = RamanSpectrum('/home/chris/Documents/DataWeiss/150424/150424_01.txt')
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150424/150424_04.txt')
    c  = RamanSpectrum('/home/chris/Documents/DataWeiss/150424/150424_10.txt')
    d  = RamanSpectrum('/home/chris/Documents/DataWeiss/150430/150430_01.txt')
    e = RamanSpectrum('/home/chris/Documents/DataWeiss/150430/150430_03.txt')
    
    ax1 = gca()
    
        
    a = removespikes(a)
    a.smooth()
    a.smooth(window_len=21, window = 'SG')
    a.autobaseline((68, 322),order = 4)
    a.autobaseline((322, 767),order = 0,join='start')
    a.autobaseline((767, 838),order = 0,join='start')
    a.autobaseline((838, 1405),order = 2,join='start')
    a.autobaseline((1405,1466),order = 0,join='start')
    a.autobaseline((1466, 1974),order = 2,join='start')
        
    b = removespikes(b)
    b.smooth(window_len=21, window = 'SG')
    b.autobaseline((68, 322),order = 4)
    b.autobaseline((322, 767),order = 0,join='start')
    b.autobaseline((767, 838),order = 0,join='start')
    b.autobaseline((838, 1405),order = 2,join='start')
    b.autobaseline((1405,1466),order = 0,join='start')
    b.autobaseline((1466, 1974),order = 2,join='start')   
#    
    c = removespikes(c)
    c.smooth(window_len=21, window = 'SG')
    c.autobaseline((68, 322),order = 4)
    c.autobaseline((322, 767),order = 0,join='start')
    c.autobaseline((767, 838),order = 0,join='start')
    c.autobaseline((838, 1405),order = 2,join='start')
    c.autobaseline((1405,1466),order = 0,join='start')
    c.autobaseline((1466, 1974),order = 2,join='start')
    
    d = removespikes(d)
    d.smooth( window_len=21, window = 'SG')
    d.autobaseline((68, 322),order = 4)
    d.autobaseline((322, 767),order = 0,join='start')
    d.autobaseline((767, 838),order = 0,join='start')
    d.autobaseline((838, 1405),order = 2,join='start')
    d.autobaseline((1405,1466),order = 0,join='start')
    d.autobaseline((1466, 1974),order = 2,join='start')
    
    e = removespikes(e)
    e.smooth(window_len=21, window = 'SG')
    e.autobaseline((68, 322),order = 4)
    e.autobaseline((322, 767),order = 0,join='start')
    e.autobaseline((767, 838),order = 0,join='start')
    e.autobaseline((838, 1405),order = 2,join='start')
    e.autobaseline((1405,1466),order = 0,join='start')
    e.autobaseline((1466, 1974),order = 2,join='start')

    

    
    
    a[:]+=200#ax1.lines[0].set_ydata(ax1.lines[0].get_ydata()+200)
    b[:]+=500
    c[:]+=800
    d[:]+=1100
    e[:]+=1500
    
    a.plot(marker = 'o',markersize = 1)
    b.plot(marker = 'o',markersize = 1)
    c.plot(marker = 'o',markersize = 1)
    d.plot(marker = 'o',markersize = 1)
    e.plot(marker = 'o',markersize = 1)
    
    ylim(0,1800)
    xlim(68,1980)
    a_list = list()
    b_list = list()
    c_list = list()
    d_list = list()
    e_list = list()
    for w in (541,628,742,1574):
        z = fitspectrum(a,(w-30,w+30),'OneGaussian', [300, w,50, 0,50])
        if z ==-1:
            print 'fit awry'
            
    
        else:
            a_list.append(z[0])
            plot(z[1], z[2])
            
    z = fitspectrum(a,(1045,1130),'TwoGaussian', [225, 1067,20,250, 1095,20, -1,50])
    if z ==-1:
        print 'fit awry'
    else:
        a_list.append(z[0])
        plot(z[1], z[2])  
    
    ###############################################################    
    for w in (623,639,793,1086,1598):
        z = fitspectrum(b,(w-30,w+30),'OneGaussian', [300, w,50, 0,850])
        if z ==-1:
            print 'fit awry'
        else:
            b_list.append(z[0])
            plot(z[1], z[2])

    ###############################################################    
    for w in (144,206,499,623,636,795,1088,1178,1280,1591):
        z= fitspectrum(c,(w-30,w+30),'OneGaussian', [300, w,50, 0,1150])
        if z ==-1:
            print 'fit awry'
        else:
            c_list.append(z[0])
            plot(z[1], z[2])

            
    ###############################################################    
    for w in (88,122,262,496,538,628,724,1066,1087,1177,1561):
        z= fitspectrum(d,(w-30,w+30),'OneGaussian', [300, w,50, 0,1550])
        if z ==-1:
            print 'fit awry'
        else:
            d_list.append(z[0])
            plot(z[1], z[2])

            
    ###############################################################    
    for w in (242,631,813,1083,1157,1588):
        z= fitspectrum(e,(w-30,w+30),'OneGaussian', [300, w,50, 0,50])
        if z ==-1:
            print 'fit awry'
        else:
            e_list.append(z[0])
            plot(z[1], z[2])
   
    f = open('/home/chris/Documents/DataWeiss/150430/thiophenolsonDots.txt', 'w')
    f.write('\n\nchloro\n\n')
    for i in a_list:
        f.write(str(i[0])+'\n')
    f.write('\n\ncmethyl\n\n')
    for i in b_list:
        f.write(str(i[0])+'\n')
    f.write('\n\nmethoxy\n\n')
    for i in c_list:
        f.write(str(i[0])+'\n')
    f.write('\n\ncbromo\n\n')
    for i in d_list:
        f.write(str(i[0])+'\n')
    f.write('\n\ncfluoro\n\n')
    for i in e_list:
        f.write(str(i[0])+'\n')
    f.close()

#    fitspectrum(spectrum, rnge, func, guess)
    return (a_list,b_list,c_list,d_list,e_list)
    
def thiophenolreffits():
    e =  RamanSpectrum('/home/chris/Documents/DataWeiss/150430/150430_08.txt')
    d = RamanSpectrum('/home/chris/Documents/DataWeiss/150430/150430_14.txt')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150424/150424_06.txt')
    c = RamanSpectrum('/home/chris/Documents/DataWeiss/141014/4_methoxythiophenol.spe')
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/141007/1_methylbenzenethiol.spe')   
        
    
    
    a.normalize()
    
    b.normalize()
    
    c.normalize()
    
    d.normalize()
    e.normalize()
    
    ax1 = gca()

    a[:]+=0
    b[:]+=1.1
    c[:]+=2.2
    d[:]+=3.3
    e[:]+=4.4
    
    a.plot(marker = 'o',markersize = 1)
    b.plot(marker = 'o',markersize = 1)
    c.plot(marker = 'o',markersize = 1)
    d.plot(marker = 'o',markersize = 1)
    e.plot(marker = 'o',markersize = 1)
    
    ylim(0,6)
    xlim(68,1980)
    a_list = list()
    b_list = list()
    c_list = list()
    d_list = list()
    e_list = list()
    for w in (541,628,734,1180,1300):
        print w
        z = fitspectrum(a,(w-30,w+30),'OneGaussian', [0.5, w,50, 0,50])
        if z ==-1:
            print 'fit awry'
        else:
            a_list.append(z[0])
            plot(z[1], z[2])
    z = fitspectrum(a,(1540,1600),'TwoGaussian', [1,1, 1561,1574,10,10,0,0])
    if z ==-1:
        print 'fit awry'
    else:
        a_list.append(z[0])
        plot(z[1], z[2])         
    z = fitspectrum(a,(1020,1120),'FourGaussian', [1,1,1,1,1060,1074, 1090, 1100,20,20,20,20,0,0])
    if z ==-1:
        print 'fit awry'
    else:
        a_list.append(z[0])
        plot(z[1], z[2])  
        
    
    ###############################################################    
    for w in (647,803,1104,1191,1220,1598):
        print w
        z = fitspectrum(b,(w-30,w+30),'OneGaussian', [0.5, w,50, 0,1.1])
        if z ==-1:
            print 'fit awry'
        else:
            b_list.append(z[0])
            plot(z[1], z[2])
   
    ###############################################################    
    for w in (206,499,623,636,795,1088,1178,1280,1591):
        print w
        print type(c)
        z= fitspectrum(c,(w-30,w+30),'OneGaussian', [1, w,50, 0,2.2])
        if z ==-1:
            print 'fit awry'
        else:
            c_list.append(z[0])
            plot(z[1], z[2])
   
            
    ###############################################################    
    for w in (262,496,538,628,724,1066,1087,1177,1561):
    
        print w
        z= fitspectrum(d,(w-30,w+30),'OneGaussian', [1, w,50, 0,3.3])
        if z ==-1:
            print 'fit awry'
        else:
            d_list.append(z[0])
            plot(z[1], z[2])
  
    ###############################################################    
    for w in (242,631,813,1083,1157,1588):
        print w
        z= fitspectrum(e,(w-30,w+30),'OneGaussian', [1, w,50, 0,4.4])
        if z ==-1:
            print 'fit awry'
        else:
            e_list.append(z[0])
            plot(z[1], z[2])
            
     ###################################################       
    f = open('/home/chris/Documents/DataWeiss/150430/thiophenolsRefs.txt', 'w')
    f.write('\n\nchloro\n\n')
    for i in a_list:
        f.write(str(i[0])+'\n')
    f.write('\n\ncmethyl\n\n')
    for i in b_list:
        f.write(str(i[0])+'\n')
    f.write('\n\nmethoxy\n\n')
    for i in c_list:
        f.write(str(i[0])+'\n')
    f.write('\n\ncbromo\n\n')
    for i in d_list:
        f.write(str(i[0])+'\n')
    f.write('\n\ncfluoro\n\n')
    for i in e_list:
        f.write(str(i[0])+'\n')
    f.close()
    
#    fitspectrum(spectrum, rnge, func, guess)
    return (a_list,b_list,c_list,d_list,e_list)
    

def May7():
    figure()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150507/150507_01.txt')
    a.normalize()
    a.plot()
    ics('/home/chris/Orca/CdTP_bridge/CdTP_bridgeDFT.out',normalize = True)
    title('thiophenol')
    
    figure()    
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150507/150507_03.txt')
    a.normalize()
    a.plot()
    i = RamanSpectrum('/home/chris/Documents/DataWeiss/150508/150508_02.txt')
    i[:]/=1200
    i.smooth()
    i.autobaseline((70,450),leaveout=(70,340),order = 4)
    i.autobaseline((450,1650),order = 2, join='start')
    i.plot()
    ics('/home/chris/Orca/CdClTP/CdClTP.out',normalize = True,labelpeaks = False)
    
    figure()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150507/150507_06.txt') ## bromocomplex
    a.autobaseline((70,450),leaveout=(70,340),order = 4)
    a.autobaseline((450,1650),order = 2, join='start')
    a.normalize()
    a.plot()
    i = RamanSpectrum('/home/chris/Documents/DataWeiss/150508/150508_08.txt')  ## bromo on dots
    i[:]/=1200
    i.smooth()
    i.autobaseline((70,450),leaveout=(70,340),order = 4)
    i.autobaseline((450,1650),order = 2, join='start')
    i.plot()
    ics('/home/chris/Orca/CdBrTP/CdBrTP.out',normalize = True,labelpeaks = False)
    
    
    
    return 0