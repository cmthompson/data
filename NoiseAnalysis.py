# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:00:46 2015

@author: chris
"""
from RamanTools3 import *

def n():  ########### remove interference noise from 4_light
   
    width = 100
    r = RamanSpectrum('/home/chris/Documents/DataWeiss/150304/4_light.SPE')
    
    lambdaa = 10**7/(10**7/514.5-array(r.index))
    r.index = pandas.Float64Index(lambdaa)
    r=r.iloc[90:900]
    
    
    for i in range(5):
        r = smooth(r)
    
    fit = polyfit(array(r.index), r.values, 4)
    r/=polyeval(fit,array(r.index))
    
    
    xs = array(r.index)/1000
    
    def funct(x,A,m,b,phase):return (1+(A)*cos((x+m*x**2)*pi/b+phase))**2
    guess = [0.005,1,0.005,0.1]  
    
    plot(xs,funct(xs,*guess),'k')
    plot(xs,r.values)
    x = scipy.optimize.curve_fit(funct,xs,r.values,guess)
    print x[0]
    plot(xs,funct(xs,*x[0]))
    
    
   
    
    return 0
    
def n2(name):  ########### remove interference noise from 3_light
   
   
    r = RamanSpectrum(name)
    
    lambdaa = 10**7/(10**7/514.5-array(r.index))
    r.index = pandas.Float64Index(lambdaa)
    r=r.iloc[300:800]
  
    weights = zeros(r.size)
    
    for i in range(5):
        r = smooth(r)
    

    fit = polyfit(array(r.index), r.values, 5)
    r/=polyeval(fit,array(r.index))
   
    
    
    
    xs = array(r.index)/1000
    
    plot(xs, r.values)

 
    def funct(p):return sum((r.values-0.0017*cos(700/(xs+4*xs**2)+p)-1)**2)
    bnds = ((0,2*pi),)   
    pguess = float(scipy.optimize.minimize(funct,pi,bounds = bnds)['x'])
    
    pguess=0
    #plot(xs,0.0017*cos(600/(xs+1*xs**2) +pguess)+1 ,'r')
    
    
    
    #def funct(x,A,m,b):return (A/1000.0)*cos(d*100/(a*x**2+b*x**2+c)+pguess)+1   #def funct(x,A,b,phase):return (1+(A)*cos(x*pi/b+phase))**2
    #guess = [1.7,4,6.0]#
    
    def funct(x,A,a,b,c):return (A/1000.0)*cos(1000/(a*x**2+b*x+c))+1   #def funct(x,A,b,phase):return (1+(A)*cos(x*pi/b+phase))**2
    guess = [1.7,3,1,0]#
    
    

   # plot(xs,funct(xs, *guess),'k')
    
    x = scipy.optimize.curve_fit(funct,xs,r.values,guess)
    print x[0]
    plot(xs,funct(xs,*x[0]),label=os.path.basename(name),color = gca().lines[-1].get_color())
    print sum((r.values-funct(xs,*x[0]))**2)
    
    legend()
   
    return 0
    
def nall():
    clf()
    names  = ['/home/chris/Documents/DataWeiss/150304/2_ light.SPE',
              '/home/chris/Documents/DataWeiss/150304/3_ light.SPE','/home/chris/Documents/DataWeiss/150304/4_light.SPE',
              '/home/chris/Documents/DataWeiss/150304/7_1800 grating.SPE',
              '/home/chris/Documents/DataWeiss/150304/8.SPE']
    a = pandas.Series([],[])
    for name in names:
        print name
        n2(name)
        
    
    
    
#    xs = array(r.index)/1000
#    
#    Aguess = max(r)
#
#    print Aguess
#    def funct(p):return sum((r.values- (Aguess)*cos(72*10/(xs)+p) - 1)**2)
#    bnds = ((0,2*pi))
#    p =float(scipy.optimize.minimize(funct,(2.2,),bounds=bnds).x)
#    print 'phase', p
#    
#    clf()
#    
#    def funct(x,A,a,b,c):return (A/1000)*cos(a*x**2+b*x+c)+1   #def funct(x,A,b,phase):return (1+(A)*cos(x*pi/b+phase))**2
#    guess = [Aguess*1000,-416,502,-146]  #guess = [0.0017,0.002,p]#
#   
#    x = scipy.optimize.curve_fit(funct,xs,r.values,guess)
#    print 'for concatenated series'
#    print x[0], p
#    plot(xs,funct(xs,*x[0]),'r',label=os.path.basename(name))
#    plot(xs,r.values,'r.')
#    
#    legend()
    
    return 0

def ncon():
    a = RamanSpectrum('/home/chris/Documents/DataWeiss/150304/2_ light.SPE')
    b=RamanSpectrum('/home/chris/Documents/DataWeiss/150304/3_ light.SPE')
    c = RamanSpectrum('/home/chris/Documents/DataWeiss/150304/4_light.SPE')
    d = pandas.concat([a,b,c])
    d=d.sort_index()
    d.plot()
    


    return d
def fitsome():
    x = array([532.1,
        534.6,
        536.8,
        539.6,
        542.5,
        545.8,
        548.7,
        552,
        555.5,
        559,
        562.7,
        566.8,
        570.7,
        574.7,
        578.9,
        583,
        587.3])
    x/=1000
    y = array([2.5,
        2.5,
        2.2,
        2.8,
        2.9,
        3.3,
        2.9,
        3.3,
        3.5,
        3.5,
        3.7,
        4.1,
        3.9,
        4,
        4.2,
        4.1,
        4.3,
        ])
    
    fit = polyfit(x,y,2)
    print fit
    return 0


    
def m():
    
    r = RamanSpectrum('/home/chris/Documents/DataWeiss/150304/2_ light.SPE')
    v = smooth(r)
    for i in range(5):
        v=smooth(v)
    v.plot(marker='s')
    r.plot()
    return 0
def o():
    figure()

    names  = ['/home/chris/Documents/DataWeiss/150304/2_ light',
          '/home/chris/Documents/DataWeiss/150304/3_ light',
              '/home/chris/Documents/DataWeiss/150304/4_light',
              '/home/chris/Documents/DataWeiss/150304/5_647',
              '/home/chris/Documents/DataWeiss/150304/6_647',
              '/home/chris/Documents/DataWeiss/150304/7_1800 grating',
              '/home/chris/Documents/DataWeiss/150304/8']
              
    ax1=subplot(221)
    ax2=subplot(222)
    ax3 = subplot(223)
    
    finenoisearray = ndarray((0,1024))
    for name in names:
        a = loadtxt(name+'.txt',delimiter=',', )[:,2:]
        a-=130
        
        freq = array(RamanSpectrum(name+'.SPE').index)
        ax3.plot(freq,a[0])
        if '647' in name:
            lambdaa = 10**7/(10**7/647.1-freq)
        else:
            lambdaa = 10**7/(10**7/514.5-freq)
            m = array([mean(a,axis=0)])
            m/=mean(m)
            
            finenoisearray = append(finenoisearray,m,axis=0)
       
        ax1.plot(lambdaa,a[0])
        
        ax2.plot(a[0])
    ax2.errorbar(range(1024), mean(finenoisearray,axis = 0), yerr=std(finenoisearray,axis=0)/2)
    ax1.legend(list(os.path.basename(r) for r in names))
    return 0
   
   

def darkcounts():  ### shows the digitized nature of dark spectra and very low signal spectra due to ADC gain setting.  
    figure()
    a = loadtxt('/home/chris/Documents/DataWeiss/150304/dark measurement 100 scans 1s each.txt',delimiter=',', unpack=True)
    b = loadtxt('/home/chris/Documents/DataWeiss/150304/18_ dark specrtrum 10s accum.txt',delimiter=',', unpack=True)
    c = loadtxt('/home/chris/Documents/DataWeiss/150304/19_ dark specrtrum 50s accum.txt',delimiter=',', unpack=True)
    d = loadtxt('/home/chris/Documents/DataWeiss/150304/20_ dark specrtrum 100s accum.txt',delimiter=',', unpack=True)
    counts = array(list((mean(i.flatten()) for i in [a,b,c,d])),'s')
    print counts
    plot([1,10,50,100],counts)
    title('total dark baseline for increasing accumulation times')
    xlabel('accumulation time (s)')
    ylabel('dark baseline (counts)')

    
    a = loadtxt('/home/chris/Documents/DataWeiss/150304/10_dark2.txt',delimiter=',', unpack=True)
    b = loadtxt('/home/chris/Documents/DataWeiss/150304/11_dark gain on 3.txt',delimiter=',', unpack=True)
    c = loadtxt('/home/chris/Documents/DataWeiss/150304/15_10s accums gain 1.txt',delimiter=',', unpack=True)
    d = loadtxt('/home/chris/Documents/DataWeiss/150304/16 10s accums gain 3.txt',delimiter=',', unpack=True)
    print list((i.shape for i in [a,b,c,d]))
    
    figure()
    av = mean(a,axis=1)
    subplot(421)
    bins = arange(-.105,.115,0.0025)
    hist(a[100]/mean(a[100])-1,bins=bins )
   
    print 'standard deviation 1 s acuumulations, gain=1:',std(a[100]/mean(a[100])-1)
    xlim(-.105,0.105)
    subplot(422) 
    hist(a[150]/mean(a[150])-1,bins=bins )
    print 'standard deviation 1 s acuumulations, gain=1:',std(a[150]/mean(a[150])-1)
    xlim(-.105,0.105)
    
    
    
    av = mean(b,axis=1)
    subplot(423)
    
    hist(b[100]/mean(b[100])-1,bins=bins )
    xlim(-.105,0.105)
    print 'standard deviation 1 s acuumulations, gain=3:',std(b[100]/mean(b[100])-1)
    subplot(424) 
    hist(b[150]/mean(b[150])-1,bins=bins )
    print 'standard deviation 1 s acuumulations, gain=3:',std(b[150]/mean(b[150])-1)
    xlim(-.105,0.105)
    
    av = mean(c,axis=1)
    subplot(425)
    
    hist(c[100]/mean(c[100])-1,bins=bins )
    xlim(-.105,0.105)
    print 'standard deviation 10 s acuumulations, gain=3:',std(c[100]/mean(c[100])-1)
    subplot(426) 
    hist(c[150]/mean(c[150])-1,bins=bins )
    print 'standard deviation 10 s acuumulations, gain=1:',std(c[150]/mean(c[150])-1)
    xlim(-.105,0.105)
    
    av = mean(b,axis=1)
    subplot(427)
    
    hist(d[100]/mean(d[100])-1,bins=bins )
    xlim(-.105,0.105)
    print 'standard deviation 10 s acuumulations, gain=3:',std(d[100]/mean(d[100])-1)
    subplot(428) 
    hist(d[150]/mean(d[150])-1,bins=bins )
    print 'standard deviation 10 s acuumulations, gain=3:', std(d[150]/mean(d[150])-1)
    xlim(-.105,0.105)
    
    return 0
    
    
    
    
    
def something(l):
    theta = 3.14159/4  #angle between grating and CCD
    a = 1E6/1100  ## size of grating groove
    beta = 3.14159/4 ## angle of incidence on the grating
    d = 30000    # thickness of layer
    
    print 'central wavelength = ', 
    #return cos(d/(l*(cos(theta)*sqrt(1-((l-a*sin(beta))/a)**2)+sin(theta)*(l-a*sin(beta))/a)))
    return cos(d/(l*cos(theta-arcsin((l-a*sin(beta))/a))))
    
def grating():
    theta = 3.14159/4  #angle between grating and CCD
    a = 1E6/1100  ## size of grating groove
    beta = 3.14159/4 ## angle of incidence on the grating
    d = 30000    # thickness of layer
    l = arange(200,600,0.01)
    plot(l,something(l))
    plot(l,cos(theta-arcsin((l-a*sin(beta))/a)))
    return 0
    
    
def v2(name):  ########### remove interference noise from 3_light
   
   
    r = RamanSpectrum(name)
    
    
    r=r.iloc[300:1000]
  
    weights = zeros(r.size)
    
    for i in range(5):
        r = smooth(r)
    
    fit = polyfit(array(r.index), r.values, 5)
    r/=polyeval(fit,array(r.index))
    
    
    
    xs = array(r.index)
    
    Aguess = max(r)

    print Aguess
    try:
        pguess = arccos(r[0]/Aguess)
    except:
        pguess = pi
    def funct(p):return sum((r.values- (Aguess)*cos(xs*pi/70+p)-1)**2)
    bnds = ((0,2*pi))
    pguess =float(scipy.optimize.minimize(funct,(pguess,),bounds=bnds).x)
    print 'phase', pguess
    
    
    
    def funct(x,b,A):return (A/1000)*cos(x*pi/b+pguess)+1   #def funct(x,A,b,phase):return (1+(A)*cos(x*pi/b+phase))**2
    guess = [70,max(r)*1000]  #guess = [0.0017,0.002,p]#
   
    x = scipy.optimize.curve_fit(funct,xs,r.values,guess)
    print x[0]
    plot(xs,funct(xs,*x[0]),label=os.path.basename(name))
    plot(xs,r.values,'.',color = gca().lines[-1].get_color())
    
    legend()
    return 0
    
def vall():
    clf()
    names  = ['/home/chris/Documents/DataWeiss/150304/3_ light.SPE',
              '/home/chris/Documents/DataWeiss/150304/4_light.SPE',
              '/home/chris/Documents/DataWeiss/150304/7_1800 grating.SPE',
              '/home/chris/Documents/DataWeiss/150304/8.SPE']
    for name in names:
        print name
        v2(name)
    return 0
    
def takeout(a,centers = (92,452,479,595,663,941),demo = False):
    
    w =11
    def function(x,A1,w1,G1,b): return b - A1*exp(-(x-w1)**2/G1)
    for center in centers:
        
       
        x = array(a.index[center-w:center+w])
        y = a.values[center-w:center+w]
        
        guess= [a.values[center]-a.values[center+w],float(a.index[center]),20,float(a.values[center+w])]
      
        try:
            result = scipy.optimize.curve_fit(function,x,y,guess)
        except RuntimeError:
            tkMessageBox.showerror('Fit Failed')
            return 0
        z = list(result[0])
        if demo:
       
            print z
            plot(x,function(x,*z))
        
        if abs(z[2]<100):
            a.iloc[center-7:center+7]-=function(x,*z)-z[-1]

    return 0