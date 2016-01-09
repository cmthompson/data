# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 11:05:14 2015

@author: chris
"""
from ramanTools.RamanSpectrum import *

def Mar24():  ########### Raman of older red dots.  These have polystyrene or toluene on them. 
    figure()
    j = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_22.txt')
    k = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_23.txt')
    l = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_26.txt')
    m = add_RamanSpectra(j,k)
    m=add_RamanSpectra(m,l)
    m=autobaseline(m,(0,3300),order = 0)
    m=smooth(m)
    m.plot(label='NativeLigands')
    
    j = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_24.txt')
    k = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_25.txt')
    l = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_26.txt')
    m = add_RamanSpectra(j,k)
    m=add_RamanSpectra(m,l)
    m=autobaseline(m,(0,3300),order = 0)
    
    m=smooth(m)
    m.plot(label='NativeLigands')
    
    j = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_28.txt')
    k = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_29.txt')
    
    m = add_RamanSpectra(j,k)
    m=autobaseline(m,(0,3300),order = 0)
    m=smooth(m)
    
    m.plot(label='NativeLigands')
    
    ref = RamanSpectrum('/home/chris/Documents/DataWeiss/141007/Liquid sample corrected-spectrum of toluene.txt')
    ref.plot(label = 'toluene')
    legend()
    title('Native Ligand')
    
    
    figure()
    j = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_34.txt')
    m=smooth(j)
    m.plot(label='MeOTP')
    
    
    j = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_35.txt')
    k = RamanSpectrum('/home/chris/Documents/DataWeiss/150324/NativeLigand_Red_36.txt')
   
    m = add_RamanSpectra(j,k)
   
    m=smooth(m)
    m.plot(label='MeOTP')
    
    title('MeOTP treated')
    
    
    return 0

def Mar26():
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150326/orangedot-nativeligand_20.txt')
    print type(a)
    a.smooth()
    
    a.autobaseline((400,520),order =0)
    a.autobaseline((520,1756),order = 4)
    a.values[:]*=10
    a.plot(label = '633 nm')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150326/orangedot-nativeligand_21.txt')
    a = smooth(a)
    a = autobaseline(a,(2482,3600),4)
    a*=10
    a.plot(label = '633 nm')
    b = CdODPARef-2597
    b.plot(label = 'reference')
    
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150326/orangedot-nativeligand_10.txt')
    a = smooth(a)
    a-=161
    a*=25
   
    a.plot(label = '785 nm')
    
    return 0
    
def etch():
    from RamanTools3 import RamanSpectrum
    os.chdir('/home/chris/Documents/DataWeiss/150326/etchegoin')
    data = zeros((0,1024))
    frequency = zeros(data.shape)
    
    l = os.listdir(os.curdir)

    l.sort()
    for x in l:
        if 'notes' not in x and '.txt' in x:
            a = loadtxt(x,unpack = True)
            data = append(data,array([a[1]]), axis = 0)
            frequency = append(frequency, array([a[0]]),axis = 0)
    print data.shape[0], 'spectra averaged'
    plot(frequency[:,512], 's')
    figure()
    ##########now the data is in the proper form.
    
    tup = copy(data)
   
    pad = zeros((tup.shape[0],512))
    tup = append(pad,tup,axis = 1)
    tup = append(tup,pad,axis = 1)   #### pad with zeroes
    Snoise = array([mean(tup,axis = 0)]*tup.shape[0])
   
    if False:
        plot(tup[0])
        plot(Snoise[0])
        plot(tup[1],'s')
        
        plot(tup[0])
        return 0
    tup2 = copy(tup)
    for i in range(tup.shape[0]):
        
        z = mean(tup[i])/mean(Snoise[0])
        
        tup[i]-=z*Snoise[0]
        tup[i] = roll(tup[i],i)  ###### correct the pixel offset
        tup2[i] = roll(tup2[i],i)  # rolling tup2 but not subtracting noise
       
   
    
    
    
    tup_av = mean(tup[:,512:1536],axis = 0)
    tup2_av = mean(tup2[:,512:1536],axis = 0)
    tup2_av-=8800
    
    tup_av = SGsmooth(arange(1024),tup_av, width = 11, order = 3)
    
    
    
    plot(arange(512,1536),tup_av,label='etchegoin') 
    plot(arange(512,1536),tup2_av,label='averaged') 
    plot(Snoise[0]-8900,label ='ccd bias')
    return 0
    
def Mar27():
    os.chdir('/home/chris/Documents/DataWeiss/150327/150327/etchegoin1')
    data = zeros((0,1024))
    frequency = zeros(data.shape)
    l = os.listdir(os.curdir)
    l.sort()
    for x in l:
        if 'notes' not in x and '.txt' in x:
            a = RamanSpectrum(x)
            a = removespikes(a)
            data = append(data,array([a.values]), axis = 0)
            frequency = append(frequency, array([a.index]),axis = 0)
    print data.shape[0], 'spectra averaged'   
    tup = copy(data)
    pad = zeros((tup.shape[0],512))
    tup = append(pad,tup,axis = 1)
    tup = append(tup,pad,axis = 1)   #### pad with zeroes
    tup2 = tup
    for i in range(tup.shape[0]):
        tup2[i] = roll(tup2[i],i)  # rolling tup2 but not subtracting noise
   
    tup2_av = mean(tup2[:,512:1536],axis = 0)
    spectrum = RamanSpectrum(pandas.Series(tup2_av, frequency[0]))
    spectrum = autobaseline(spectrum,(763,1620), order = 4)
    spectrum=smooth(spectrum)
    spectrum.plot(color='k')
    #################################################
    os.chdir('/home/chris/Documents/DataWeiss/150327/150327/etchegoin2')
    data = zeros((0,1024))
    frequency = zeros(data.shape)
    l = os.listdir(os.curdir)
    l.sort()
    for x in l:
        if 'notes' not in x and '.txt' in x:
            a = RamanSpectrum(x)
            try:
                a = removespikes(a)
            except:
                print a
                return 0
                
            data = append(data,array([a.values]), axis = 0)
            frequency = append(frequency, array([a.index]),axis = 0)
    print data.shape[0], 'spectra averaged'   
    tup = copy(data)
    pad = zeros((tup.shape[0],512))
    tup = append(pad,tup,axis = 1)
    tup = append(tup,pad,axis = 1)   #### pad with zeroes
    tup2 = tup
    for i in range(tup.shape[0]):
        tup2[i] = roll(tup2[i],i)  # rolling tup2 but not subtracting noise
   
    tup2_av = mean(tup2[:,512:1536],axis = 0)
    spectrum = RamanSpectrum(pandas.Series(tup2_av, frequency[0]))
    spectrum = autobaseline(spectrum,(124,900), order = 4)
    spectrum=smooth(spectrum)
    spectrum.plot(color = 'k')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150327/150327/MeOTP treated new _73.txt')
    a=autobaseline(a,(1200,2000),order = 4)
    a = smooth(a)
    
    
    a.plot(color='k')
    ##################################################
    

   
  
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150327/red dots new sample/orangedot-Nativeligand_26.txt')
   
    a=autobaseline(a,(0,1700),order=4)
  
    a=smooth(a)
    a/=100
    a+=100
    a.plot(color = 'r')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150327/red dots new sample/orangedot-Nativeligand_24.txt')
    a=autobaseline(a,(75,1700),order = 4)
    a = smooth(a)
    a/=10
    a+=100
    a.plot(color='r')

   
    ylim(-50,250)
    xlim(100,1800)
    annotate('native ligand', (10000,250),fontsize=30,color = 'r')
    annotate('methoxy treated', (1000,200), fontsize = 30,textcoords = 'axes fraction',color = 'k')
    return 0
    
def temp():
    os.chdir('/home/chris/Documents/DataWeiss/150327/150327/etchegoin1')
    data = zeros((0,1024))
    frequency = zeros(data.shape)
    l = os.listdir(os.curdir)
    l.sort()
    for x in l:
        if 'notes' not in x and '.txt' in x:
            a = RamanSpectrum(x)
            a = removespikes(a)
            data = append(data,array([a.values]), axis = 0)
            frequency = append(frequency, array([a.index]),axis = 0)
    print data.shape[0], 'spectra averaged'   
    tup = copy(data)
    pad = zeros((tup.shape[0],512))
    tup = append(pad,tup,axis = 1)
    tup = append(tup,pad,axis = 1)   #### pad with zeroes
    tup2 = tup
    for i in range(tup.shape[0]):
        tup2[i] = roll(tup2[i],i)  # rolling tup2 but not subtracting noise
   
    tup2_av = mean(tup2[:,512:1536],axis = 0)
    spectrum = RamanSpectrum(pandas.Series(tup2_av, frequency[0]))
   # spectrum = autobaseline(spectrum,(763,1620), order = 4)
   # spectrum=smooth(spectrum)
    spectrum.plot(color='k')
    #################################################
    os.chdir('/home/chris/Documents/DataWeiss/150327/150327/etchegoin2')
    data = zeros((0,1024))
    frequency = zeros(data.shape)
    l = os.listdir(os.curdir)
    l.sort()
    for x in l:
        if 'notes' not in x and '.txt' in x:
            a = RamanSpectrum(x)
            try:
                a = removespikes(a)
            except:
                print a
                return 0
                
            data = append(data,array([a.values]), axis = 0)
            frequency = append(frequency, array([a.index]),axis = 0)
    print data.shape[0], 'spectra averaged'   
    tup = copy(data)
    pad = zeros((tup.shape[0],512))
    tup = append(pad,tup,axis = 1)
    tup = append(tup,pad,axis = 1)   #### pad with zeroes
    tup2 = tup
    for i in range(tup.shape[0]):
        tup2[i] = roll(tup2[i],i)  # rolling tup2 but not subtracting noise
   
    tup2_av = mean(tup2[:,512:1536],axis = 0)
    spectrum = RamanSpectrum(pandas.Series(tup2_av, frequency[0]))
    
    spectrum.plot(color = 'k')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150327/150327/MeOTP treated new _73.txt')
    
    
    
    a.plot(color='k')
    return 0

def Mar31():
    
    subplot(221)
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150331/150331_02.txt')
    a = autobaseline(a,(141,1700),order = 4)
    a = smooth(a)+50
    #b = RamanSpectrum('/home/chris/Documents/DataWeiss/150331/150331_03.txt')
    #c =add_RamanSpectra(a,b)
    a.plot(color = 'k')
    

    
    ylim(0,200)
    xlim(140,1700)
    
    subplot(222)
    
    a.plot(color = 'k')
    xlim(2700,3100)
    ylim(50,100)
    
    subplot(223)
    (RamanSpectrum('/home/chris/Documents/DataWeiss/140918/9_CdMeOTP.SPE')/10-1800).plot(color = 'b')
    CdODPARef.plot(color = 'r')
    xlim(100,1700)
    ylim(0,5000)
    
    subplot(224)
    (RamanSpectrum('/home/chris/Documents/DataWeiss/140918/7_CdMeOTP.SPE')-40000).plot(color='b')
    (CdODPARef+10000).plot(color = 'r')
    xlim(2700,3100)

def Apr6():
    subplot(221)
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150406/150406_02.txt')
    a.plot(color = 'k')
    
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150406/150406_03.txt')
    a.plot(color = 'k')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150406/150406_04.txt')
    a.plot(color = 'k')
    
    ylim(1000,3500)
    xlim(100,1700)
    
    subplot(222)
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150406/150406_01.txt')
    a.plot(color = 'k')
    xlim(2700,3100)
    
    subplot(223)
    (RamanSpectrum('/home/chris/Documents/DataWeiss/140918/9_CdMeOTP.SPE')/10-1800).plot(color = 'b')
    CdODPARef.plot(color = 'r')
    xlim(100,1700)
    ylim(0,5000)
    
    subplot(224)
    (RamanSpectrum('/home/chris/Documents/DataWeiss/140918/7_CdMeOTP.SPE')-40000).plot(color='b')
    (CdODPARef+10000).plot(color = 'r')
    xlim(2700,3100)

    return 0
    
def Apr7():
    a = loadtxt('/home/chris/Documents/DataWeiss/150407/150407edit.csv',delimiter = ',', unpack=True,skiprows=1)
    peakmax= argmin(abs(571-a[0]))
    plot(a[0], a[2])
    figure()
    for i in range(1,9):
        a[i]-=min(a[i])
        a[i]/=a[i,peakmax]
        plot(a[0],a[i],label=str(i))
        
#    plot(a[0],a[1]-min(a[1]),label='1')
#    plot(a[0],a[2]-min(a[2]),label='2')
#    plot(a[0],a[3]-min(a[3]),label='3')
#    plot(a[0],a[4]-min(a[4]),label='4')
#    plot(a[0],a[5]-min(a[5]),label='5')
#    plot(a[0],a[6]-min(a[6]),label='6')
#    plot(a[0],a[7]-min(a[7]),label='7')
#    plot(a[0],a[8]-min(a[8]),label='8')
    legend()
    vlines(571,0,1)
    return 0
    
def Apr8UVVis():
    a = loadtxt('/home/chris/Documents/DataWeiss/150408/UVVis/150408_edited.csv',delimiter = ',', unpack=True,skiprows=1)
    peakmax= argmin(abs(571-a[0]))

    for i in range(3,7):
        a[i]-=min(a[i])
        a[i]/=a[i,peakmax]
    fig = gcf()
    ax1 = fig.add_subplot(121)
    
    ax1.plot(a[0],a[3],label='cleaned 4x',color = 'k',linewidth=4)
    ax1.plot(a[0],a[5],label='cleaned4x + 178eq MeOTP 45 min',color = 'r',linewidth=4)
    legend()
    ax1.set_ylim(-0.1,2)
    ax1.set_xlim(450,800)
    ax1.set_title('570 nm dots cleaned 4x + 178eq MeOTP')
    ax1.set_xlabel('wavelength (nm)')
    ax1.set_ylabel('absorbance (a.u.)')
    
    ax3 = fig.add_axes((0.3,0.6,0.15,0.2))
    ax3.plot(a[0],a[3],label='cleaned 4x',color = 'k',linewidth=4)
    ax3.plot(a[0],a[5],label='cleaned4x + 178eq MeOTP 45 min',color = 'r',linewidth=4)
    ax3.set_ylim(0.8,1.1)
    ax3.set_xlim(530,610)
    
    
    ax2 = fig.add_subplot(122)
    ax2.plot(a[0],a[4],label='cleaned 5x',color = 'k',linewidth=4)
    ax2.plot(a[0],a[6],label='cleaned5x + 178eq MeOTP 45 min',color = 'r',linewidth=4)
    ax2.set_ylim(-0.1,2)
    ax2.set_xlim(450,800)
    ax2.legend()
    ax2.set_title('570 nm dots cleaned 5x + 178eq MeOTP')
    ax2.set_xlabel('wavelength (nm)')
    ax2.set_ylabel('absorbance (a.u.)')
    ax4 = fig.add_axes((0.73,0.6,0.15,0.2))
    ax4.plot(a[0],a[4],label='cleaned 4x',color = 'k',linewidth=4)
    ax4.plot(a[0],a[6],label='cleaned4x + 178eq MeOTP 45 min',color = 'r',linewidth=4)
    ax4.set_ylim(0.8,1.1)
    ax4.set_xlim(530,610)
    
    
    return 0
    
def Apr8Raman():
    os.chdir('/home/chris/Documents/DataWeiss/150408')
    fig = figure()
    a = RamanSpectrum('150408_15.txt')
    a = autobaseline(a, (200,1700),leaveout=(200,300),order=4)
    a+=800
    
    b = RamanSpectrum('150408_02.txt')
    b = autobaseline(b,(200,1700),leaveout=(200,300), order = 4)
    
    (normalize(MeOTPRef,(0,10000))*4000+2000).plot(color ='b',linewidth=2)
    a.plot(color = 'k',linewidth = 2)
    b.plot(color = 'r', linewidth = 2)
    
    ylim(-500,10000)
    xlim(200,1675)
    
    
    
    legend(['MeOTP ref', 'MeOTP treated','Native ligand only'])
    
    ylabel('Raman Intensity (a.u.)')
    xlabel('Raman Shift (cm$^{-1}$)')
    
    figure()
    title('Washing')
    a = RamanSpectrum('150408_11.txt')
    b = RamanSpectrum('150408_02.txt')
    a = autobaseline(a, (200,1700),leaveout=(200,300),order=4)
    b = autobaseline(b,(200,1700),leaveout=(200,300), order = 4)
   
    (normalize(ODPARef,(0,10000))*4000+2000).plot(color ='b',linewidth=2, label='ODPA Ref')
    a.plot(color = 'r',label='washed 5x')
    b.plot(color = 'k',label='washed 4x')
    a= fitspectrum(b,(900,1150),'SixGaussian', [200,200,200,200,200,200,950,990,1026,1064,1087,1115,10,10,10,10,10,10,1,-100])
    plot(a[1],a[2],linewidth =3,label='Cdenriched fit')
    legend()
    ylabel('Raman Intensity (a.u.)')
    xlabel('Raman Shift (cm$^{-1}$)')
    return 0
    
def Apr8Raman_forVictor():
    os.chdir('/home/chris/Documents/DataWeiss/150408')
    fig = figure(figsize=(12,6))
    subplot(121)
    a = RamanSpectrum('150408_15.txt')
    a = autobaseline(a, (200,1700),leaveout=(200,300),order=4)
    a+=800
    
    
    
    b = RamanSpectrum('150408_02.txt')
    b = autobaseline(b,(200,1700),leaveout=(200,300), order = 4)
    
    (normalize(MeOTPRef,(0,10000))*4000+2000).plot(color ='b',linewidth=2,label = 'thiophenolate reference')
    a.plot(color = 'k',linewidth = 2)
    b.plot(color = 'r', linewidth = 2)
    
    ylim(-500,10000)
    xlim(740,1675)
 
    arrowprops={'width':1,'headwidth':3,'color':'k'}
    ylabel('Raman Intensity (a.u.)')
    xlabel('Raman Shift (cm$^{-1}$)')
    annotate('C-S-H bend', (913,2830),xytext = (913,3300), xycoords = 'data',arrowprops = arrowprops,horizontalalignment='center' )
    subplot(122)
    e = RamanSpectrum('150408_13.txt')
    e.autobaseline((2500,3600),leaveout=(200,300),order=2)
    e.autobaseline((2500,2800),leaveout=(200,300), order = 1,join='end')
    e+=800
    
    
    
    f = RamanSpectrum('150408_03.txt')
    f.autobaseline((2500,3600),leaveout=(200,300), order = 2)
    f.autobaseline((2500,2800),leaveout=(200,300), order = 1,join='end')
    
    (normalize(MeOTPRef,(0,10000))*4000+2000).plot(color ='b',linewidth=2)
    e.plot(color = 'k',linewidth = 2, )
    f.plot(color = 'r', linewidth = 2)
    annotate('S-H stretch', (2560,3370),xytext = (2600,4500), xycoords = 'data',arrowprops = arrowprops,horizontalalignment='center' )
   
    ylim(-500,10000)
    xlim(2500,3200)
    legend(['thiophenol reference','CdSe thiophenolate-treated','CdSe native ligand only'])
    
    ylabel('Raman Intensity (a.u.)')
    xlabel('Raman Shift (cm$^{-1}$)')

   # savetxt()

    return 0
    
def Apr15():
    def gauss(x,A,G,m,b):return A*exp(-(1090-x)**2/G)+m*x+b
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150415/150415_15.txt')
    a = SPIDcorrect(a)
    noise = calc_noise(a,(900,1000))
    start = argmin(abs(1065-array(a.index)))
    end = argmin(abs(1115-array(a.index)))
    x = array(a.index[start:end])
    y = a.values[start:end]
    guess = [10,15,-1,y[0]+1000]
    
    peak = scipy.optimize.curve_fit(gauss, x, y, guess)
    print peak[0]
    
    
    a.plot(color = 'k')
    #plot(x,gauss(x,*peak[0]))
    #
    print 'signal to noise =', sqrt(peak[0][0]/noise)
    ylim(27000,37000)
    ylabel('Raman Intensity (a.u.)')
    xlabel('Raman Shift (cm$^{-1}$)')
    
    
    ax3 = gcf().add_axes((0.6,0.6,0.25,0.25))
    a.plot(color='k', ax = ax3)
    ax3.annotate('S/N: '+str(1.66), (1080, 28700), textcoords = 'data', size = 18)
    ax3.set_ylim(27500,29300)
    ax3.set_xlim(900,1200)
    return 0
   
    
def Apr20():
    os.chdir('/home/chris/Documents/DataWeiss/150420')
    a = add_RamanSpectra(RamanSpectrum('150420_03.txt'), RamanSpectrum('150420_04.txt'))
    
    a =SPIDcorrect633(a)
    a = smooth(a)
    a.plot()    
    a = autobaseline(a, (700,1200), order = 3)
    a = autobaseline(a, (100,698), order=4)
    #a = autobaseline(a, (300,500), order=2)
    #a = autobaseline(a, (100,300), order=2, leaveout=(190,220))

    
    figure()
    
    i_list = ['150420_04.txt',
              '150420_05.txt',
              '150420_06.txt',
              '150420_07.txt',
              '150420_09.txt',
              '150420_10.txt',
              '150420_11.txt',
              '150420_12.txt',
              '150420_13.txt',
              '150420_14.txt',
              '150420_15.txt',
              '150420_16.txt']
    multlist = [1,0.30, 0.400, 0.100, 0.600, 0,0.20, 0.400, 0.10, 1.00,1.00,1.00]
    ax1 = gcf().add_subplot(121)
    ax2 = gcf().add_subplot(122)
    a= RamanSpectrum(pandas.Series(zeros((2048,)), linspace(50,1600, 2048)))
    def gauss(x,A,m,b):return A*exp(-(1093-x)**2/10)+m*x+b
    ston_list = list()
    
    for i in i_list:
        b = RamanSpectrum(i)
        b.plot(ax= ax1)
        #b = removespikes(b)
        b = SPIDcorrect633(b)
        b[:]*=multlist[i_list.index(i)]
        s = calc_noise(b, (800,900))
        start = argmin(abs(1070-array(b.index)))
        end = argmin(abs(1110-array(b.index)))
        x = array(b.index[start:end])
        y = b.values[start:end]
        guess = [10,-1,y[0]+1000]
        
        peak = scipy.optimize.curve_fit(gauss, x, y, guess)
        ax1.plot(x,gauss(x,*peak[0]))
        
       
        print i, peak[0][0], s
        
        
       
        a = add_RamanSpectra(a, b)
        
        s = calc_noise(a, (800,900))
        start = argmin(abs(1070-array(a.index)))
        end = argmin(abs(1110-array(a.index)))
        x = array(a.index[start:end])
        y = a.values[start:end]
        guess = [10,-1,y[0]+1000]
        
        peak = scipy.optimize.curve_fit(gauss, x, y, guess)
        print i, peak[0][0], s
        ston_list.append(peak[0][0]/s)
        a.plot(ax = ax2)
        ax2.plot(x,gauss(x,*peak[0]))
    ax1.legend(i_list)    
        
    a = smooth(a,window_len=7, window = 'SG')
  
    a.plot()
    figure()
    plot(ston_list)
        
    return 0
    
def Apr20UVVis():
    a = loadtxt('/home/chris/Documents/DataWeiss/150420/150420.csv',delimiter = ',', unpack=True,skiprows=1)
    peakmax= argmin(abs(571-a[0]))
    
    fig = gcf()
    ax1=fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    for i in range(1,8):
        a[i]-=min(a[i])
        a[i]/=a[i,peakmax]
        if i == 7: 
            lstyle = 'k--'
        else:
            lstyle='-'
        ax1.plot(a[0],a[i],lstyle, linewidth=3)
        ax2.plot(a[0],a[i],lstyle, linewidth=3)
   
    
    ax1.legend(['0 eq', '50 eq','100 eq', '200 eq', '500 eq', '1000 eq', '500 eq octanethiol'])
    ax1.set_ylim(0,2)
    ax1.set_xlim(450,800)
    ax1.set_title('570 nm dots titrated with MeOTP')
    ax1.set_xlabel('wavelength (nm)')
    ax1.set_ylabel('absorbance (a.u.)')
    ax1.vlines(633, 0,2,linewidth = 3)
    
    ax2.set_ylim(0.97,1.02)
    ax2.set_xlim(562,583)
    ax2.set_title('570 nm dots titrated with MeOTP')
    ax2.set_xlabel('wavelength (nm)')
    ax2.set_ylabel('absorbance (a.u.)')

    return 0
    
def Apr22():
    os.chdir('/home/chris/Documents/DataWeiss/150422')


    
    figure()
    
    i_list = ['150422_03.txt',
              '150422_05.txt',
              '150422_06.txt',
              '150422_07.txt',
              '150422_09.txt']
#              '150422_11.txt',
#              '150422_12.txt',
#              '150422_13.txt',
#              '150422_14.txt',
#              '150422_15.txt',
#              '150422_16.txt']
    ax1 = gcf().add_subplot(121)
    ax2 = gcf().add_subplot(122)
    a= RamanSpectrum(pandas.Series(zeros((2048,)), linspace(50,1600, 2048)))
    def gauss(x,A,m,b):return A*exp(-(1093-x)**2/10)+m*x+b
    ston_list = list()
    
    for i in i_list:
        b = RamanSpectrum(i)
        b = removespikes(b)
        b = SPIDcorrect633(b)
        if i ==i_list[0]:
            b*=10
        elif i == i_list[2]:
            b.iloc[700]-=10000
        b.plot(ax= ax1)
        
        

        a = add_RamanSpectra(a, b)
        a-=min(a.values)
        
        s = calc_noise(a, (800,900))
        start = argmin(abs(1070-array(a.index)))
        end = argmin(abs(1110-array(a.index)))
        x = array(a.index[start:end])
        y = a.values[start:end]
        guess = [10,-1,y[0]+1000]
        
        peak = scipy.optimize.curve_fit(gauss, x, y, guess)
        print i, peak[0][0], s
        ston_list.append(peak[0][0]/s)
        a.plot(ax = ax2)
        ax2.plot(x,gauss(x,*peak[0]))
    ax1.legend(i_list)    
        
    a = smooth(a,window_len=7, window = 'SG')
  
    a.plot()
    figure()
    plot(ston_list)
        
    return 0
    