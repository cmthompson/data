# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/chris/.spyder2/.temp.py
"""
from ramanTools.RamanSpectrum import *
from UVVistools import *
from scipy.special import erf
from numpy.polynomial.polynomial import polyder,polyval
from numpy.polynomial.polynomial import polyfit as pf
import weissdatavariables



def findderivative(x,y,derivorder=1):

    window_len=25
#    if spectrum.ndim != 1:
#        raise ValueError, "smooth only accepts 1 dimension arrays."
#
#    if spectrum.size < window_len:
#        raise ValueError, "Input vector needs to be bigger than window size."

        
    

    retval = numpy.ndarray(y.shape)
    for i in range(y.size):
        
        i_min = max(0,i-window_len/2)
        i_max = min( i+window_len/2-1, y.size)
        fit = pf(x[i_min:i_max+1],y[i_min:i_max+1],min(5,window_len))
        retval[i] = polyval(x[i],polyder(fit,derivorder),tensor=False) 
        
    return retval
    

def Nov6TCSPC(guess):
    """TCSPC results from PPA capped CdS dots"""
   # os.chdir('/media/cybertron_box/Chris Thompson/TCSPC Data')
    os.chdir('151106/TCSPC Data')
    files = ['11-6-15-chris-cds-exchanged-pH10-375exc-100psbin---433nm_em_621p83t',
             '11-6-15-chris-cds-375exc-50psbin--run2-423nm_em_3797p62t',
             '11-6-15-scattercell-375ex-750em-213p43t', 
             '11-6-15-scattercell-375ex-750em-68p68t',
             '11-6-1--5-chris-cds-375exc-50psbin-423nm_em_1265p45t',
             '11-6-15-chris-cds-exchanged-pH10-375exc-100psbin---433nm_em_171p27t',
             '11-6-15-chris-cds-exchanged-375exc-100psbin---430nm_em_841p41t']
    

    pH10_1 = loadtxt(files[0],unpack = True,)
    oleate_2 = loadtxt(files[1],unpack = True)
    IRF_1 = loadtxt(files[2],unpack = True)
    IRF_2 = loadtxt(files[3],unpack = True)
    oleate_1 = loadtxt(files[4],unpack = True)
    pH7 = loadtxt(files[5],unpack = True)
    pH10_2 = loadtxt(files[6],unpack = True)
    
     

    
    pH10_1[1]+=pH10_2[1]
   
    
    oleate_1[1]+=oleate_2[1]
    
    
    IRF_1[1]+=IRF_2[1]
    for i  in [oleate_1, pH7, pH10_1,IRF_1]:
        zero = numpy.mean(i[1,10:argmax(i[1])-10],axis = 0)
        
        i[1]-=zero
        i[1]/=max(i[1])
        
    
#    t0x = argmax(IRF_1[1])
#    t0 = IRF_1[0,t0x]
    
        
    def gaussian(x,A,w,x0):return (A/sqrt(2*pi)/w)*exp(-(x-x0)**2/(2*w**2))
        
    fitIRF = scipy.optimize.curve_fit(gaussian, IRF_1[0], IRF_1[1],[1,1,45])[0]
    fixt0=fitIRF[2]
  
    w = fitIRF[1]
    
    
    
    
    
    
    w0 =w# w/(2*sqrt(2*log(2)))
    print w0, fixt0
        
   
    def fitfunction(x,*args):
        y = numpy.zeros(x.shape)
        t0=args[0]
        for z in range((len(args)-1)/2):
            
            A1 = args[z*2+1]
            tau1 = args[z*2+2]
           # y =(A1/2)*exp(w0^2/(2*tau1^2) - (t-t0)/tau1) * (1 - erf((w0^2-tau1*(t-t0))/(w0*tau1*sqrt(2)))); 
            y+=(A1/2)*exp(w0**2/(2*tau1**2) - (x-t0)/tau1)* (1 - erf((w0**2-tau1*(x-t0))/(w0*tau1*sqrt(2))))  # / (exp(w0**2/(2*tau1**2))* (1 - erf((w0)/(tau1*sqrt(2)))))
        return y
        

        
    

#    
    legend(['oleate','PPA pH7','PPA pH10', 'IRF'])
    
    for i  in [IRF_1]:
       # plot(i[0],fitfunction(i[0],*guess))
        fit = scipy.optimize.curve_fit(fitfunction, i[0], i[1],[44,1,0.1])[0]
        print fit
        plot(i[0],fitfunction(i[0],*fit),linewidth=3)
    figure()
    for r  in range(1,4):
        
        subplot(310+r)
        i = list([oleate_1, pH7, pH10_1])[r-1]
        print numpy.mean(numpy.diff(i[0]))
        plot(i[0],i[1],'.')
       # plot(i[0],fitfunction(i[0],*guess))
        fit = scipy.optimize.curve_fit(fitfunction, i[0], i[1],guess,maxfev=3000)[0]
        
        xlim(40,60)
        ylim(-0.05,1.05)
       
        print list(("%.2f" % i for i in fit))
        plot(i[0],fitfunction(i[0],*fit),linewidth=3)
        l = list([fit[0],fit[1],fit[2]])
        plot(i[0],fitfunction(i[0],*l),linewidth=1)
        
        l = list([fit[0],fit[3],fit[4]])
        plot(i[0],fitfunction(i[0],*l),linewidth=1)
        
        if len(fit)>6:
            l = list([fit[0],fit[5],fit[6]])
            plot(i[0],fitfunction(i[0],*l),linewidth=1)
        
        tau = array(fit[2::2])
        A = array(fit[1::2])
        tau_ave = (sum(A*tau**2)/sum(A*tau))
        print 'tau ave', tau_ave
        
        annotate(['oleate','PPA pH7','PPA pH10'][r-1],(0.8,0.8),xycoords='axes fraction')
    xlabel('time (ns)')
    
    

    
    return 0

#fixt0=0
#def fitfunction(x,*args):
#    y = numpy.zeros(x.shape)
#    t0=fixt0
#    for z in range((len(args))/2):
#        print A1,tau1
#        A1 = args[z*2+0]
#        tau1 = args[z*2+1]
#       # y =(A1/2)*exp(w0^2/(2*tau1^2) - (t-t0)/tau1) * (1 - erf((w0^2-tau1*(t-t0))/(w0*tau1*sqrt(2)))); 
#        y+=(A1/2)*exp(w0**2/(2*tau1**2) - (x-t0)/tau1)* (1 - erf((w0**2-tau1*(x-t0))/(w0*tau1*sqrt(2))))
#    return y
    

def Nov18():
    a =loadtxt('/home/chris/Dropbox/DataWeiss/151118/Dots with untreated alumina.csv', unpack = True, delimiter = ',', skiprows = 1,usecols = (0,1,7,9,11))
    s = pf(a[0][0:150],a[4][0:150],1)
    a[4]-=polyval(a[0],s)
    for i in a[array([1,3,4])]:#[1:]:
        i-=i[150]
        r = findpeak(a[0],i,(405,420))
        print r
        
        i-=i[150]
        i/=r[1]
        plot(a[0], i)
    
    legend(['as prep', 'alumina 4 mg', 'alumina 5.8 mg'])
    ylabel('absorbance')
    xlabel('wavelength')
    return 0
    
def Nov23():
    figure()
    a = loadtxt('/home/chris/Dropbox/DataWeiss/151110/ExchangeYieldforPPAexchange.csv', unpack = True, delimiter = ',', skiprows = 2,usecols = (0,1))
    for i in a[1:]:
        i[:]-=i[0]
        r = findpeak(a[0],i,(405,420))
       
        x = argmin(abs(a[1,0:190]-r[1]/2))
        print 'FWHM', a[0,x]-r[0]
        print r
        print CdSconc(r),'M'
        plot(a[0],i)
    print '---------------'
    figure()
    a =loadtxt('/home/chris/Dropbox/DataWeiss/151123/151123PPAdots_forICP.csv', unpack = True, delimiter = ',', skiprows = 1,usecols = (0,1,2,3,17,18,19,21,22,23))
    x = argmin(abs(a[0]-450))    
    for i in a[1:]:
        i[:]-=i[x]
        r = findpeak(a[0],i,(405,420))
        print r
       # print CdSconc(r)*5,'M'
        plot(a[0],i)
    print 'subractiing from 450'
   
    print 'using fixed epsilon---------------'
    eps =  476708.893203
#    for i in a[1:]:
#        i[:]-=i[x]
#        r = findpeak(a[0],i,(405,420))
#        conc = r[1]/eps/0.2
#        #print conc,'M'
#        plot(a[0],i)
    return 0
    
    
def Nov25UVVis():
   
    figure()
    a = loadtxt('/home/chris/Dropbox/DataWeiss/151125/151125ppatitratin.csv', unpack = True, delimiter = ',', skiprows = 2)
    for i in a[1:]:
        x = argmin(abs(a[0]-450))
        i[:]-=i[x]
        r = findpeak(a[0],i,(405,420))
        i[:]/=r[1]
        print r
        plot(a[0],i)
    with open ('/home/chris/Dropbox/DataWeiss/151125/151125ppatitratin.csv','rb') as f:
        r = f.readline().split(',')[1:]
        
        f.close()
        
    legend(r)
    return 0

 