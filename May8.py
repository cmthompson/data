# -*- coding: utf-8 -*-
"""
Created on Sun May 10 20:32:31 2015

@author: chris
"""
from ramanTools.RamanSpectrum import *
def May8():
    os.chdir('/home/chris/Documents/DataWeiss/150508')
    x = RamanSpectrum('150508_01.txt')
    cl = RamanSpectrum('150508_02.txt')
    cl.values[:]*=2
    h =  RamanSpectrum('150508_04.txt')
    c8 = RamanSpectrum('150508_05.txt')
    ome = RamanSpectrum('150508_07.txt')
    br = RamanSpectrum('150508_08.txt')
    py = RamanSpectrum('150508_09.txt')
    brwash = RamanSpectrum('150508_10.txt')
    me = RamanSpectrum('150508_11.txt')
    clf()
    ax1 = subplot(121)
    for i in (cl,h,ome,br,me):
        i.smooth()
        i.autobaseline((70,450),leaveout=(70,340),order = 4)
        i.autobaseline((450,1650),order = 2, join='start')
        i.plot()
    legend([ 'cl', 'h', 'ome','br','me'])
    ax2 = subplot(122)
    for i in (x,py,c8):
        i.smooth()
        i.autobaseline((50,400),leaveout=(70,340),order = 4)
        i.autobaseline((400,1650),order = 3, join='start')
        i.plot()
    legend(['x','py','c8'])
    
    
    
    
 
    return 0
    
def May14():
    clf()
    s = RamanSpectrum('/home/chris/Documents/DataWeiss/150514/150514_12.txt')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150514/150514_13.txt')
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150514/150514_14.txt')
    c = add_RamanSpectra(a,b)
   # c = add_RamanSpectra(s,c)
   
    d = RamanSpectrum('/home/chris/Documents/DataWeiss/150514/150514_15.txt')
    v = RamanSpectrum('/home/chris/Documents/DataWeiss/150514/150514_16.txt')
    z=0.15
    
    e = subtract_RamanSpectra(c, d*z)
    
    #e.smooth()
    l = subtract_RamanSpectra(c,v*z)
    #l.smooth()
    e = RamanSpectrum(e.append(l))
    e = autobaseline(e,(180,277),order = 4)
    e = autobaseline(e,(277,1700),order = 4,join='start')
    #j = subtract_RamanSpectra(c,RamanSpectrum(pandas.Series(f[2],f[1])))
    e.plot()
    
    #d.plot()
    #v.plot()
    
    
   # CdMeOTPRef.index = array(CdMeOTPRef.index)-5
    (CdMeOTPRef/120).plot()
    (MeOTPRef/240).plot()
  
  
    return 0
    
def May15():
    clf()
    #s = RamanSpectrum('/home/chris/Documents/DataWeiss/150515/150515_02.txt')
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150515/150515_03.txt')
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150515/150515_04.txt')
    c = add_RamanSpectra(a,b)
    
    c = autobaseline(c,(300,1700),order = 4)
    c.smooth()
    c.plot(label='Cd-enriched',color = 'k')
    a= fitspectrum(c,(900,1150),'SixGaussian', [200,200,200,200,200,200,950,990,1026,1064,1087,1115,10,10,10,10,10,10,1,-100])
    plot(a[1],a[2],linewidth =3,color = 'k')
    
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150515/150515_07.txt')
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150515/150515_08.txt')
    c = add_RamanSpectra(a,b)
    
    c = autobaseline(c,(300,1700),order = 4)
    c.smooth()
    c.plot(label='stoichiometric',color = 'r')
    
    
   # CdMeOTPRef.index = array(CdMeOTPRef.index)-5
   # (CdMeOTPRef/120).plot()
   # (MeOTPRef/240).plot()
    (RamanSpectrum('/home/chris/Documents/DataWeiss/150511/150511_01.txt')/4).plot()
    (CdODPARef/10-100).plot()
   # ics('/home/chris/Orca/Successful/CdMeOTP/CdMeOTP.out')
   # ics('/home/chris/Orca/CdTP_bridge/CdTP_bridgeDFT.out',color='r')
    a= fitspectrum(c,(900,1150),'SixGaussian', [200,200,200,200,200,200,950,990,1026,1064,1087,1115,10,10,10,10,10,10,1,-100])
    plot(a[1],a[2],linewidth =3,color = 'r')
    return 0
    
def May16():
    
    
    c = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_08.txt')
    c.values[:]*=3
    c = autobaseline(c,(300,1700),order =6)
    c.smooth()
    c.plot(label='Cd-enriched')
    a= fitspectrum(c,(900,1150),'SixGaussian', [200,200,200,200,200,200,950,990,1026,1064,1087,1115,10,10,10,10,10,10,1,-100])
    plot(a[1],a[2],linewidth =3,label='Cdenriched fit')
    
#    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_08.txt')
#    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_07.txt')
#    c = add_RamanSpectra(a,b)
#    
#    c = autobaseline(c,(300,1700),order = 4)
#    c.smooth()
#    c.plot(label='stoichiometric')
#    
    
   # CdMeOTPRef.index = array(CdMeOTPRef.index)-5
   # (CdMeOTPRef/120).plot()
   # (MeOTPRef/240).plot()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_01.txt')
    b = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_02.txt')
    c = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_03.txt')*4
    d = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_05.txt')
    e = RamanSpectrum('/home/chris/Documents/DataWeiss/150516/150516_06.txt')
    a = add_RamanSpectra(a,b)
   
    a = add_RamanSpectra(a,c)

    
    
    a = add_RamanSpectra(a,d)
    
    a = add_RamanSpectra(a,e)
    
    a.values[:]/=10
    a.plot(label='pieces')
  
   # ics('/home/chris/Orca/Successful/CdMeOTP/CdMeOTP.out')
   # ics('/home/chris/Orca/CdTP_bridge/CdTP_bridgeDFT.out',color='r')
    #a= fitspectrum(a,(900,1150),'SixGaussian', [200,200,200,200,200,200,950,990,1026,1064,1087,1115,10,10,10,10,10,10,1,-100])
   # plot(a[1],a[2],linewidth =3,label= 'piecesfit')
    return 0
    
def Voigt(x,A,w0,G,g):
    return A*exp(-(x-w0)**2/G)/(3.14159*((x-w0)**2+g**2))

def May18():
    r = RamanSpectrum('/home/chris/Documents/DataWeiss/150518/150518_03b.txt')
    r.autobaseline((100,763),order =3)
    r.autobaseline((763,836),order =1,join='start')
    r.autobaseline((836,1600),order =3,join='start')
    r.plot()
    
    r = RamanSpectrum('/home/chris/Documents/DataWeiss/150518/150518_05b.txt')
    r.autobaseline((100,763),order =3)
    r.autobaseline((763,836),order =1,join='start')
    r.autobaseline((836,1650),order =3,join='start')
    r.plot()
    return 0

def today():
    a = '/home/chris/Documents/DataWeiss/150518/'
    for i in range(1,9):
        r = RamanSpectrum(a+'150517_0'+str(i)+'.txt')
        v = 10**7/array(r.index)-10**7/785
        r = RamanSpectrum(pandas.Series(r.values,-1*v))
        r.to_csv(a+'150518_0'+str(i)+'b.txt')

    return 0    
    

def May20():
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150520/150520_02.txt')
#    def xGaussian(x,*guess):
#        numpeaks = (len(guess)-2)/3
#        y = guess[-2]*x/1000+guess[-1]
#        for i in range(numpeaks):
#            y+= guess[i]*exp(-(x-guess[i+numpeaks])**2/guess[i+2*numpeaks])
#        return y
    n_guess = [100,100,100,100,100,100,100,100,720,790,785,811,820,848,865,885,10,10,10,10,10,10,10,10,0,600]
    
    a.plot()
    b = fitspectrum(a,(700,908),'xGaussian',n_guess)
    plot(b[1],b[2])
    return 0
    
def May21():
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150521/150521stoic_dots.CSV')
    a[:]-=0.2
    n_guess = [0.05,0.05,0.05,0.1,0.05,0.05,0.05,0.05,0.05,0.05,
               1015, 1048,1075,1100,1159,1169,1184,1204,1221,1246,
               20,20,20,40,20,40,20,20,20,20,
               0,0.0]
    
    a.plot(color = 'k')
    print a.nearest(1100)
    print a.nearest(1300)
    b = fitspectrum(a,(1000,1260),'xGaussian',n_guess)
    print b.params[0]
    plot(b.x,b.y,'r')
    print len(b.x)
    print len(b.peaks)
    #plot(b.x,b.peaks[0])
    for p in b.peaks:
        plot(b.x,p,'b')
        #pass

    b = ((RamanSpectrum('/home/chris/Documents/DataWeiss/150520/150520_02.txt')-500)/10000)   
    b.plot()
    return 0

def May21b():
    clf()
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150520/150520_02.txt')
    a[:]-=0.2
    n_guess = [0.05,0.05,0.05,0.05,0.05,0.05,0.05,
               1068,1100,1169,1184,1204,1221,1246,
               20,20,20,20,20,20,20,
               0,0.21]
    
    a.plot(color = 'k')
    CdODPARef.autobaseline((200,1700),order = 1)
    CdODPARef.plot()
    
   # b = RamanSpectrum('/home/chris/Documents/DataWeiss/150408/150408_02.txt')
   # b = autobaseline(b,(200,1700),leaveout=(200,300), order = 4)
    #b.plot()
    return 0
    
def May27():
    clf()
    ax1 = subplot(121) 
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150527/150527CdOPA2_from_May8_heated.CSV')
#    p = [892,915,929,945,962,970,989,1005,1030,1046,1038,1061,1072,1095,1108,1135,1122]
#    n_guess = [0.5]*len(p)+[0.2] + p+[1075] + [25]*len(p)+[10000] +[0,0]
#               
#    
#    a.plot(color = 'k',marker = '.')
#    
#    b = fitspectrum(a,(860,1150),'xGaussian',n_guess)
#    print b.params[0]
#    plot(b.x,b.y,'r')
# 
#    
#    for p in b.peaks:
#        plot(b.x,p,'b')
#        #pass)
    p = [915,929,945,955,964,968,989,1005,1030]
    n_guess = [0.5]*len(p) + p + [25]*len(p) +[0,0]
               
    
    a.plot(color = 'k',marker = '.')
    
    b = fitspectrum(a,(923,1011),'xGaussian',n_guess)
    print b.params[0]
    plot(b.x,b.y,'r')
 
    
    for p in b.peaks:
        plot(b.x,p,'b')
        #pass)
    
    ax2 = subplot(122)   
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150527/150527CdOPA_ala_Cao_heated.CSV')
    p = [890,916,955,965]#,983,990,1008,1030,1049,1061,1071,1090,1098,1110,1120]
    A = [0.2,0.05,0.4,0.1]#,0.2,0.2,0.1,0.1,0.4,0.12,0.3,0.3,0.4,0.2,0.2]    
    n_guess = A + p + [25]*len(p) +[0,0.03]
    
    a.plot(color = 'k',marker = '.')
    
    b = fitspectrum(a,(860,930),'xGaussian',n_guess)
    print b.params[0]
    plot(b.x,b.y,'r')
 
    
    for p in b.peaks:
        plot(b.x,p,'b')
        #pass)
        
    p = [965,983,990,1008,1030]#,1049,1061,1071,1090,1098,1110,1120]
    A = [0.1,0.2,0.2,0.1,0.1]#0.4,0.12,0.3,0.3,0.4,0.2,0.2]    
    n_guess = A + p + [25]*len(p) +[0,0.03]
    a.plot(color = 'k',marker = '.')
    b = fitspectrum(a,(860,1017),'xGaussian',n_guess)
    print b.params[0]
    plot(b.x,b.y,'r')
   
    for p in b.peaks:
        plot(b.x,p,'b')
        #pass)
    p = [1049,1061,1071,1090,1098,1110,1120]
    A = [0.4,0.12,0.3,0.3,0.4,0.2,0.2]    
    n_guess = A + p + [25]*len(p) +[0,0.03]
    a.plot(color = 'k',marker = '.')
    b = fitspectrum(a,(1017,1150),'xGaussian',n_guess)
    print b.params[0]
    plot(b.x,b.y,'r')
   
    for p in b.peaks:
        plot(b.x,p,'b')
        #pass)
        
        
    return 0