# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 14:07:03 2015

@author: chris
"""
from ramanTools.RamanSpectrum import *
def spinning():  ## test of optical damage decreased by spinning

    os.chdir('/home/chris/Documents/DataWeiss/150209')
    figure()
    a =  ['1_sample not moving.SPE', 
        '2_sample not moving.SPE',
        '3_sample not moving 5 min irradiation.SPE',
        '4_sample not moving 5 min irradiation.SPE', 
        '5_spinning.SPE',
        '6_spinning.SPE', 
        '7_spinning 2 min.SPE',
        '8_spinning 3 min.SPE', 
        '9_spinning 4 min.SPE',
        '10_spinning 5 min.SPE',
        '11_spinning 10 min.SPE', '12_stopped spinning 0 min.SPE', '13_ 1min later.SPE', '14_ 2min later.SPE', '15_ 3min later.SPE', '16_ 4min later.SPE', '17_long scan.SPE', '18_red dots.SPE', '19.SPE']
    areas = list()    
    for i in a[0:17]:
        
        r = RamanSpectrum(i)
        
        print r.name
        areas.append(r.calc_area((180,220)))
        
    plot([0,0.5,5,5.5],areas[0:4]/areas[0],'s-')
    plot([0,1,2,3,4,5,10],areas[4:11]/areas[4],'s-')
    plot([0,1,2,3,4],areas[11:16]/areas[11],'s-')
    legend(['not spun','while spinning','stopped spinning'])
    return areas
    
    
def N2():  ## test of optical damage decreased by spinning

    os.chdir('/home/chris/Documents/DataWeiss/150210')
    figure()
    print os.listdir('.')
    
    a = [ '1_under N2 0 min.SPE', '1_under N2 1 min.SPE', '3_under N2 2 min.SPE', 
    '4_under N2 3 min.SPE', '5_under N2 4 min.SPE', '6_under N2 5 min.SPE', 
    '7_under N2 6 min.SPE', '8_under N2 7 min.SPE', '9_under N2 8min.SPE',
    '10_under N2 9min.SPE', '11_under N2 10min.SPE', '12_ spinning under N2 0 min.SPE', 
    '13_ spinning under N2 1min.SPE', '14_ spinning under N2 2min.SPE', '15_ spinning under N2 3min.SPE',
    '16_ spinning under N2 4min.SPE', '17_ spinning under N2 5min.SPE', 
    '18_ spinning under N2 6min.SPE', '19_ spinning under N2 7min.SPE',
    '20_ spinning under N2 8min.SPE', '21_ spinning under N2 9min.SPE',
    '22 spinning under N2 10min.SPE', '23_after dark time.SPE',
    '24_after dark time 1min.SPE', '25_ spinning in air 0 min.SPE', 
    '26_ spinning in air 1 min.SPE', '27_ spinning in air 2 min.SPE',
    '28_ spinning in air3 min.SPE', '29_ spinning in air 4min.SPE', 
    '30_ spinning in air 5min.SPE', '31_ spinning in air 6 min.SPE',
    '32_ spinning in air 7 min.SPE', '33_ spinning in air 8 min.SPE', 
     '34_ spinning in air 9 min.SPE','34_ spinning in air 10min.SPE',
    '35_after 10 min dark.SPE', '36_after 10 min dark plus 1min.SPE',
    '37.SPE', '38_ 30 SEC.SPE', '39_ 60s.SPE',
    '40_90 s.SPE', '41_120s.SPE', '42_150s.SPE']
 
    areas = list()    
    for i in a:
        
        r = RamanSpectrum(i)
        
       
        areas.append(r.calc_area((180,220)))
    N2only = areas[0:11] 
    N2spin = areas[11:24]
    airspin = areas[24:37]
    airstill = areas[37:]
    
    N2only/=N2only[0]
    N2spin/=N2spin[0]
    airspin/=airspin[0]
    airstill/=airstill[0]
    plot(N2only,'s-')
    plot([0,1,2,3,4,5,6,7,8,9,10,20,21],N2spin,'s-')
    plot([0,1,2,3,4,5,6,7,8,9,10,20,21],airspin,'s-')
    plot(arange(0,3,0.5),airstill,'s-')
    fill_between((10,20),(0,0),(1,1),color = 'y')
    annotate('dark', (15,0.4),horizontalalignment = 'center', fontsize=24)
    legend(['n2 only','n2 + spinning','air + spinning','air only'])
    
    
    xlabel('Time (min)')
    ylabel('Phonon Mode Intensity (a.u.)')
    return areas
    
def Feb10():
    figure()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150210/43_ long scan.SPE')
    
    
    adding  = pandas.Series([NaN]*len(arange(300,1500,0.5)),arange(300,1500,0.5))
    
    d=a.append(adding)
    
    d = d.interpolate(method='index')

    d = d[arange(300,1500,0.5)]
    
    e = FourierFilter(d,width = 1100)
    
    e.plot()
    
    
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150210/44.SPE')
    c=b+a
    a.autobaseline((300,1600),order = 4)
    #a = smooth(a,window_len=9)
    
    b.autobaseline((300,1600),order = 4)
    #b=smooth(b,window_len=9)
    
    #a.plot()
    #b.plot()
    
    c= autobaseline(c,(300,1600),order = 4)
    c = smooth(c, window_len=9)
    #c.plot()
    legend(['a','d','e'])
    
   
    
    return d
    
def LB_XPS():   ##### Stephanie made some LB films of oleate capped particles

    pandas.options.display.mpl_style = False
    spot1 = loadtxt('/home/chris/Documents/DataWeiss/150211/Spot1.csv',
            delimiter = ',',
            skiprows = 1,
            unpack = True)
    spot2 = loadtxt('/home/chris/Documents/DataWeiss/150211/Spot2.csv',
            delimiter = ',',
            skiprows = 1,
            unpack = True)
            
    
            
    print 'first the spot1 then spot2'
    print '-------------------------------------'
    for i in range(2):
        if i == 0:
            spec = spot1
            title = 'spot1'
        elif i == 1:
            spec = spot2
            title = 'spot2'
       
            
        
        
        figure().suptitle(title,size = 24)
        subplot(331)
        a = RamanSpectrum(pandas.Series(spec[13],spec[12]))
        a.plot()
        cadmium_area = calc_area(a,(408,401))+calc_area(a,(416,409))
        cadmium5_2_area = calc_area(a,(408,401))
        
        start = argmin(abs(array(a.index)-408))
        end = argmin(abs(array(a.index)-401))
        print start,end
        xs =a.index[start:end]
        ys= a.values[start:end]
        
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = slope*(xs-xs[0])+ys[0]
        fill_between(xs,baseline,ys)
        annotate('Cd 3d', (0.1,0.1),xycoords = 'axes fraction')
        
        xlim((390,420))
        
        
        
        subplot(332)
        a = RamanSpectrum(pandas.Series(spec[11],spec[10]))
        silver_area = calc_area(a,(378,365))#+calc_area(a,(376,372))
        
        a.plot()
        start = argmin(abs(array(a.index)-378))
        end = argmin(abs(array(a.index)-365))
        print start,end
        xs =a.index[start:end]
        ys= a.values[start:end]
        
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = slope*(xs-xs[0])+ys[0]
        fill_between(xs,baseline,ys)
        axis('on')
        xlim((350,390))
        annotate('Ag 3d', (0.1,0.1),xycoords = 'axes fraction')
        #################
        subplot(333)
        a = RamanSpectrum(pandas.Series(spec[3],spec[2]))
        
        a.plot()
        if i == 1:
            selenium_area = calc_area(a,(56,52))
            start = argmin(abs(array(a.index)-56))
            end = argmin(abs(array(a.index)-52))
            print start,end
            xs =a.index[start:end]
            ys= a.values[start:end]
            if len(xs)==0 or len(ys)==0:
                pass
            else:
            
                slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
                baseline = slope*(xs-xs[0])+ys[0]
                fill_between(xs,baseline,ys)
        else:
            selenium_area=0
        
        xlim((45,65))
        annotate('Se 3d', (0.1,0.1),xycoords = 'axes fraction')
        #print 'selenium atoms per silver atom:', (sulfur_area/0.6666)/(silver_area/5.987)
        #### sensitivities for 54.6 degrees,  come from Handbook of XPS, p. 253
        
        subplot(334)
        a = RamanSpectrum(pandas.Series(spec[7],spec[6]))
        
        a.plot()
        xlim((155,165))
        sulfur_area = calc_area(a,(165,155))
        annotate('S 2p', (0.1,0.1),xycoords = 'axes fraction')
        
        
        subplot(335)
        a = RamanSpectrum(pandas.Series(spec[5],spec[4]))
        
        a.plot()
        start = 0#argmin(abs(array(a.index)-135))
        end = argmin(abs(array(a.index)-130))
       
        xs =a.index[start:end]
        ys= a.values[start:end]
        
        slope =(ys[-1]-ys[0])/(xs[-1]-xs[0])
        baseline = array([ys[-1]]*len(xs))
        phosphorous_area = sum((ys-baseline)[1:]*diff(xs))
        fill_between(xs,baseline,ys)
        xlim((120,140))
       # sulfur_area = calc_area(a,(175,155))
       # annotate('S 2p', (0.1,0.1),xycoords = 'axes fraction')
        
        
        subplot(336)
        a = RamanSpectrum(pandas.Series(spec[9],spec[8]))
        a.plot()
        annotate('carbon', (0.1,0.1), xycoords = 'axes fraction')
        
        
        subplot(339)
        a = RamanSpectrum(pandas.Series(spec[1],spec[0]))
        
        a.plot()
        
        
        #### sensitivities for 54.6 degrees,  come from Handbook of XPS, p. 253
     
        
        
        
        
        #print 'Cd area total', cadmium_area
        #print 'Cd 3d_5/2 only area',cadmium5_2_area
        #print 'Ag area', silver_area
        print 'Se atoms per silver atom', (selenium_area/0.853)/(silver_area/5.987)
        print 'Cd atoms per silver atom', (cadmium5_2_area/3.98)/(silver_area/5.987)
        print 'sulfur atoms per silver atom:', (sulfur_area/0.6666)/(silver_area/5.987)
        print 'P atoms per silver atom', (phosphorous_area/0.486)/(silver_area/5.987)
        print '__________________________________-'
    return 0