# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 18:34:28 2016

@author: chris
"""
from UVVistools import *

import weissdatavariables
import numpy as np

from scipy.special import erf
from ramanTools.RamanTools import modify_plot
from ramanTools.RamanSpectrum import *
from ramanTools.RamanReferences import *

def Feb5QuantumYieldsofCdSePPA():
    """Calculated quantum yields of PPA-capped CdSe dots and MPA capped CdSe in PPA buffers. 160205"""
    fig =figure()
    fig2 = figure()
    Aax1 = fig.add_subplot(221)
    ax1 = fig2.add_subplot(221)
    Aax2 = fig.add_subplot(222)
    ax2 = fig2.add_subplot(222)
    Aax3 = fig.add_subplot(223)
    ax3 = fig2.add_subplot(223)
    Aax4 = fig.add_subplot(224)
    ax4 = fig2.add_subplot(224)
    sample1 = list()
    sample2 = list()
    sample3 = list()
    sample4 = list()
    

    uvvisfile='160205/160205_CdSePPAuvvis.csv'
    fluorfolder = '160205/160205_fluor'
    rhodaminefile = '160205/160205_fluor/RhodamineB.dat'
    

    ######## Calculate Concentrations of dots transfered and the yield
    a = loadtxt(uvvisfile, skiprows = 1, usecols = (0,1,2,3,4,5,11,17,23),unpack=True,delimiter=',')
    concs=array([])
    for i in a[1:]:
        peak = findpeak(a[0],i,(525,540))
        concs = append(concs,QDconc(peak))
    concs[4:]*=(2250/250)
    print concs
    

    ## yields
    
    print "oleate dots QY:",indivCdSeQY(uvvisfile, 'ad','ae',fluorfolder +'/OleateDots2.dat', rhodaminefile,fluorescencerange = (500,615), standardfluorescencerange=(520,700),excitationwavelength=500,fluorplot=ax1,nliq=1.375)
  
    
    sample1.append(indivCdSeQY(uvvisfile, 'f','ae',fluorfolder +'/Sample1pH6.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = '',color = 'r'))
    sample1.append(indivCdSeQY(uvvisfile, 'g','ae',fluorfolder +'/Sample1pH7.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = '',color = 'r'))
    sample1.append(indivCdSeQY(uvvisfile, 'h','ae',fluorfolder +'/Sample1pH8.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = '',color = 'r'))
    sample1.append(indivCdSeQY(uvvisfile, 'i','ae',fluorfolder +'/Sample1pH9.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = '',color = 'r'))
    sample1.append(indivCdSeQY(uvvisfile, 'j','ae',fluorfolder +'/Sample1pH10.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = '',color = 'r'))
    sample1.append(indivCdSeQY(uvvisfile, 'k','ae',fluorfolder +'/Sample1pH11.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = '',color = 'r'))
    
    
    sample2.append(indivCdSeQY(uvvisfile, 'l','ae',fluorfolder +'/Sample2pH6.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = '',color = 'r'))
    sample2.append(indivCdSeQY(uvvisfile, 'm','ae',fluorfolder +'/Sample2pH7.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = '',color = 'r'))
    sample2.append(indivCdSeQY(uvvisfile, 'n','ae',fluorfolder +'/Sample2pH8.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = '',color = 'r'))
    sample2.append(indivCdSeQY(uvvisfile, 'o','ae',fluorfolder +'/Sample2pH9.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = '',color = 'r'))
    sample2.append(indivCdSeQY(uvvisfile, 'p','ae',fluorfolder +'/Sample2pH10.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = '',color = 'r'))
   #sample2.append(indivCdSeQY(uvvisfile, 'q','ae',fluorfolder +'/Sample2pH11.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'CdSe1pH11',color = 'r'))
    
    sample3.append(indivCdSeQY(uvvisfile, 'r','ae',fluorfolder +'/Sample3pH6.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = '',color = 'r'))
    sample3.append(indivCdSeQY(uvvisfile, 's','ae',fluorfolder +'/Sample3pH7.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = '',color = 'r'))
    sample3.append(indivCdSeQY(uvvisfile, 't','ae',fluorfolder +'/Sample3pH8.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = '',color = 'r'))
    sample3.append(indivCdSeQY(uvvisfile, 'u','ae',fluorfolder +'/Sample3pH9.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = '',color = 'r'))
    sample3.append(indivCdSeQY(uvvisfile, 'v','ae',fluorfolder +'/Sample3pH10.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'CdSe1pH10',color = 'r'))
    sample3.append(indivCdSeQY(uvvisfile, 'w','ae',fluorfolder +'/Sample3pH11.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'CdSe1pH11',color = 'r'))
   
    sample4.append(indivCdSeQY(uvvisfile, 'x','ae',fluorfolder +'/Sample4pH6.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH6',color = 'r'))
    sample4.append(indivCdSeQY(uvvisfile, 'y','ae',fluorfolder +'/Sample4pH7.dat', rhodaminefile,fluorescencerange = (524,590), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH7',color = 'r'))
    sample4.append(indivCdSeQY(uvvisfile, 'z','ae',fluorfolder +'/Sample4pH8.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH8',color = 'r'))
    sample4.append(indivCdSeQY(uvvisfile, 'aa','ae',fluorfolder +'/Sample4pH9.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH9',color = 'r'))
    sample4.append(indivCdSeQY(uvvisfile, 'ab','ae',fluorfolder +'/Sample4pH10.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH10',color = 'r'))
    sample4.append(indivCdSeQY(uvvisfile, 'ac','ae',fluorfolder +'/Sample4pH11.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH11',color = 'r'))
    indivCdSeQY(uvvisfile, 'ac','ae',fluorfolder +'/Sample4nobuff.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH11',color = 'r')
    
    indivCdSeQY(uvvisfile, 'ac','ae',fluorfolder +'/water.dat', rhodaminefile,fluorescencerange = (524,615), standardfluorescencerange=(520,700),excitationwavelength=500,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'CdSe1pH11',color = 'r')
    

    ax1.legend(['6','7','8','9','10','11'])
    ax2.legend(['6','7','8','9','10','11'])
    ax3.legend(['6','7','8','9','10','11'])
    ax4.legend(['6','7','8','9','10','11'])
    Aax1.legend(['6','7','8','9','10','11'])
    Aax2.legend(['6','7','8','9','10','11'])
    Aax3.legend(['6','7','8','9','10','11'])
    Aax4.legend(['6','7','8','9','10','11'])
    
    Aax1.set_ylim(0,0.2)
    Aax2.set_ylim(0,0.2)
    Aax3.set_ylim(0,0.2)
    Aax4.set_ylim(0,0.2)
    
    figure()
    pHs = [6,7,8,9,10.2,11.2]
    plot(pHs,sample1)
    plot(pHs[:-1],sample2)
    plot(pHs,sample3)
    plot(pHs,sample4)
    print 'average QYs', (array(sample1[:-1])+array(sample2[:])+array(sample3[:-1]))/3
    print (sample1[-1]+sample3[-1])/2
    print 'MPA QYs', sample4
    hlines(0.052, 6,13)
    
#    if True:######### This makes a spreadsheet for a figure for the paper
#        UVVis = loadtxt(uvvisfile,skiprows = 1, usecols = (0,10,28,29),delimiter = ',',unpack=True)
#        UVVis[1:]-=transpose([UVVis[1:,0]]*UVVis.shape[1])
#        figure()
#        for i in range(1,4):
#            a = findpeak(UVVis[0], UVVis[i],(520,540))
#            print QDconc(a)
#            UVVis[i]/=a[1]
#            plot(UVVis[0], UVVis[i])
#        
#        savetxt('../../Raul/Ligand Exchange outline/FiguresAndData/UVVisCdSe.csv',transpose(UVVis),delimiter = ',', header = 'wavelength, PPAcappedCdSepH11, MPAcappedCdSepH11, oleate capped CdSe in hexanes' )
#
#        figure()
#        b = loadtxt(fluorfolder +'/Sample1pH11.dat',usecols =(0,3), delimiter = '\t',unpack=True)
#        water = loadtxt(fluorfolder +'/water.dat',usecols =(0,3), delimiter = '\t',unpack=True)
# 
#        
#        fl = transpose([ax1.lines[-1].get_xdata(), ax1.lines[-1].get_ydata()])
#        savetxt('../../Raul/Ligand Exchange outline/FiguresAndData/FluorescenceCdSeOleate.csv',fl,delimiter = ',', header = 'wavelength, PPAcappedCdSepH11, MPAcappedCdSepH11, oleate capped CdSe in hexanes' )
#        fl = transpose([ax1.lines[-2].get_xdata(), ax1.lines[-2].get_ydata()])
#        savetxt('../../Raul/Ligand Exchange outline/FiguresAndData/FluorescenceCdSePPA.csv',fl,delimiter = ',', header = 'wavelength, PPAcappedCdSepH11, MPAcappedCdSepH11, oleate capped CdSe in hexanes' )
#        fl = transpose([ax4.lines[-2].get_xdata(), ax4.lines[-2].get_ydata()])
#        savetxt('../../Raul/Ligand Exchange outline/FiguresAndData/FluorescenceCdSeMPA.csv',fl,delimiter = ',', header = 'wavelength, PPAcappedCdSepH11, MPAcappedCdSepH11, oleate capped CdSe in hexanes' )
###        
    return (sample1, sample2, sample3, sample4)
    

    
def Feb8CdSeResonanceRaman():
    """ResonanceRaman spectra of CdSe dots with PPA ligand in water."""
    for i in range(10,17):
        a= RamanSpectrum('160208/160208_'+str(i)+'.txt')
        a.smooth()
        a.autobaseline((0,600),2)
        a.plot()
   
    modify_plot(ax=gca())
    
    return 0

    


def Feb9TCSPC(numberofexponentials):
    """TCSPC results from PPA capped CdS dots on Feb9 at pH 6 8 and 11"""
   # os.chdir('/media/cybertron_box/Chris Thompson/TCSPC Data')

    files = ['160209/2-9-16-irf-375-1mhz-scattercell',   
             '160209/2-9-16-375ex-250psbin-cds-hexanes-chris-400p855t',  
       '160209/2-9-16-375ex-250psbin-cds-aqeous+sample1-chris-1549p13t', ###  CdS1 pH 6
       '160209/2-9-16-375ex-250psbin-cds-aqeous-ph8--sample1-chris--1210p29t',  ##CdS1  pH 8
    '160209/2-9-16-375ex-250psbin-cds-aqeous-ph11--sample1-chris--599p9t',  ##CdS1  pH11
      '160209/2-9-16-375ex-250psbin-cds-aqeous-sample2-ph6-430em-chris-1617p43t', 
       '160209/2-9-16-375ex-250psbin-cds-aqeous-sample2-ph8-430em-chris-1114p9t',
    '160209/2-9-16-375ex-250psbin-cds-aqeous-sample2-ph12-430em-chris-602p585t',
    '160209/2-9-16-375ex-250psbin-cds-aqeous-ph9--sample1--430em-chris--2670p5t',
    '160209/2-9-16-375ex-250psbin-cds-aqeous-sample2-ph9-430em-chris',   ##sample2 pH9
    '160209/2-9-16-375ex-250psbin-cds-aqeous-sample1-ph10-430em-chris',
    '160209/2-9-16-375ex-250psbin-cds-aqeous-sample2-ph10-430em-chris-654p146t']
    IRF_1 = loadtxt(files[0],unpack = True)
    oleate_1 = loadtxt(files[1],unpack = True,)

    pH6 = loadtxt(files[2],unpack = True)
    pH8 = loadtxt(files[3],unpack = True)
    pH11 = loadtxt(files[4],unpack = True)
    CdS1pH9 = loadtxt(files[8],unpack = True)
    CdS1pH10 = loadtxt(files[10],unpack = True)
    CdS2pH6 = loadtxt(files[5],unpack = True)
    CdS2pH8 = loadtxt(files[6],unpack = True)
    CdS2pH11 = loadtxt(files[7],unpack = True)
    CdS2pH9 = loadtxt(files[9],unpack = True)
    CdS2pH10 = loadtxt(files[11],unpack = True)


    ######## FIT THE IRF TO A GAUSSIAN
    def gaussian(x,A,w,x0):return (A/sqrt(2*pi)/w)*exp(-(x-x0)**2/(2*w**2))
        
    fitIRF = scipy.optimize.curve_fit(gaussian, IRF_1[0], IRF_1[1],[1,1,45])[0]
    #plot(IRF_1[0], gaussian(IRF_1[0],*fitIRF))
    fixt0=fitIRF[2]
    
    w = fitIRF[1]

    w0 =w#/(2*sqrt(2*log(2)))  #### w0 is the FWHM of the gaussian used to fit the IRf
    print w0, fixt0, "IRF width and time"
 

   
    for i  in [oleate_1, pH6, pH8,pH11,CdS1pH9, CdS1pH10, CdS2pH6, CdS2pH8,CdS2pH9, CdS2pH10, CdS2pH11]:
       # i[1] = SGsmooth(i[0],i[1])
        zero = numpy.mean(i[1,10:argmax(i[1])-10],axis = 0)
     
        i[1]-=zero
        i[1]/=max(i[1])
        
 
    ## FIT THE DATA TO THE IRF CONVOLUTED WITH EXPONENTIALS
    def fitfunction(x,*args):
        """t0, A1, tau1, A2, tau2..."""
        y = numpy.zeros(x.shape)
        t0=args[0]
        for z in range((len(args)-1)/2):
            
            A1 = args[z*2+1]
            tau1 = args[z*2+2]
           
            y+=(A1/2)*exp(w0**2/(2*tau1**2) - (x-t0)/tau1)* (1 - erf((w0**2-tau1*(x-t0))/(w0*tau1*sqrt(2))))  # / (exp(w0**2/(2*tau1**2))* (1 - erf((w0)/(tau1*sqrt(2)))))
        return y

    figure()
    suptitle('sample1')
    
 
  
    errors = array([])
    timeconstants_s1 = ndarray((0,3),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        
    print '-----performing fits -----'
    for r  in range(1,7):
        
        subplot(610+r)
        if numberofexponentials==3:
            if r==0:
                guess = [44,0.05,12,0.2,6,0.04,0.5]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6,0.2,1,0.04,0.5]#,0.04, 1]
            else:
                guess = [44.3,0.05,40, 0.01,25,    0.06,   5]
        elif numberofexponentials ==2:
            if r==0:
                guess = [44,0.05,12,0.04,0.5]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6,0.04,0.5]#,0.04, 1]
            else:
                guess = [44.3,0.05,40,    0.06,   5]
        elif numberofexponentials ==1:
            if r==0:
                guess = [44,0.05,12]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6]#,0.04, 1]
            else:
                guess = [44.3,0.05,40]
        data = list([oleate_1, pH6, pH8,CdS1pH9, CdS1pH10, pH11])[r-1]
        #print numpy.mean(numpy.diff(i[0]))
        
       # plot(i[0],fitfunction(i[0],*guess))   ### plot the guess
        try:
            fit = scipy.optimize.curve_fit(fitfunction, data[0], data[1],guess,maxfev=3000)[0]
            
        except:
            fit = array(guess)
        #v=fit[1] + fit[3] +fit[5]
        xlim(-1,20)
        ylim(-0.05,1.05)
       
        print list(("%.3f" % i for i in fit))
        plot(data[0]-fit[0],data[1],'.')  #### plot the data
        plot(data[0]-fit[0],fitfunction(data[0]-fit[0],0,*fit[1:]),linewidth=3) ### plot the total fit
        errorstart = argmin(abs(42-data[0]))
        errors = np.append(errors,sqrt(sum(((data[1][errorstart:]-fitfunction(data[0][errorstart:],*fit))/data[1][errorstart:])**2)))
        
        # ratearray = array([(fit[1], fit[2]), (fit[3], fit[4]), (fit[5], fit[6])],     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray = array(list((fit[i+1], fit[i+2]) for i in arange(numberofexponentials)*2),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray=np.sort(ratearray,order = 'tau')
        for i in ratearray:
            plot(data[0]-fit[0],fitfunction(data[0]-fit[0],0,i['amp'],i['tau']),linewidth=1)
    
        

        timeconstants_s1= np.append(timeconstants_s1,ratearray)
        
        annotate(['oleate','PPA pH6','PPA pH8','pH9', 'pH10','PPA pH11'][r-1],(0.8,0.8),xycoords='axes fraction')
    timeconstants_s1=timeconstants_s1.reshape((6,-1))
    
    xlabel('time (ns)')
    

    timeconstants_s2 =ndarray((0,3),     dtype=[('amp', '<f8'), ('tau', '<f8')])
    
    figure()
    suptitle('sample2')
    for r  in range(1,7):
        
        subplot(610+r)
        if numberofexponentials==3:
            if r==0:
                guess = [44,0.05,12,0.2,6,0.04,0.5]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6,0.2,1,0.04,0.5]#,0.04, 1]
            else:
                guess = [44.3,0.05,40, 0.01,25,    0.06,   5]
        elif numberofexponentials ==2:
            if r==0:
                guess = [44,0.05,12,0.04,0.5]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6,0.04,0.5]#,0.04, 1]
            else:
                guess = [44.3,0.05,40,    0.06,   5]
        elif numberofexponentials ==1:
            if r==0:
                guess = [44,0.05,12]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6]#,0.04, 1]
            else:
                guess = [44.3,0.05,40]
        data = list([oleate_1, CdS2pH6, CdS2pH8,CdS2pH9, CdS2pH10, CdS2pH11])[r-1]
        #print numpy.mean(numpy.diff(i[0]))
        
       # plot(i[0],fitfunction(i[0],*guess))   ### plot the guess
        try:
            fit = scipy.optimize.curve_fit(fitfunction, data[0], data[1],guess,maxfev=3000)[0]
        except:
            fit = array(guess)
        #v=fit[1] + fit[3] +fit[5]
        xlim(40,60)
        ylim(-0.05,1.05)
       
        print list(("%.3f" % i for i in fit))
        plot(data[0]-fit[0],data[1],'.')  #### plot the data
        plot(data[0]-fit[0],fitfunction(data[0],0,*fit[1:]),linewidth=3) ### plot the total fit
        errorstart = argmin(abs(42-data[0]))
        errors = np.append(errors,sqrt(sum(((data[1][errorstart:]-fitfunction(data[0][errorstart:],*fit))/data[1][errorstart:])**2)))
        
        
        ratearray = array(list((fit[i+1], fit[i+2]) for i in arange(numberofexponentials)*2),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray=np.sort(ratearray,order = 'tau')
        for i in ratearray:
            plot(data[0],fitfunction(data[0],0,i['amp'],i['tau']),linewidth=1)
        

        timeconstants_s2= np.append(timeconstants_s2,ratearray)
        annotate(['oleate','PPA pH6','PPA pH8','pH9','pH10','PPA pH11'][r-1],(0.8,0.8),xycoords='axes fraction')
    timeconstants_s2=timeconstants_s2.reshape((6,-1))
    
    print 'average error', np.mean(errors)
    xlabel('time (ns)')
    return (timeconstants_s1, timeconstants_s2)



def fitplot(info):
    (timeconstants_s1, timeconstants_s2)=info
    clf()
    tfig= gcf()
    tax=tfig.add_subplot(211)
    aax=tfig.add_subplot(212)
    suptitle('timeconstants')
#    tax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s1_fast[1:],'r')
#    tax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s2_fast[1:],'b')
#    tax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s1_slow[1:],'r')
#    tax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s2_slow[1:],'b')

    for i in range(timeconstants_s1.shape[1]):
        print timeconstants_s1[1:,i]['tau']
        tax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s1[1:,i]['tau'],'rs-')
        tax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s2[1:,i]['tau'],'bs-')

        aax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s1[1:,i]['amp'],'rs-')
        aax.plot([7.0,7.8,8.4,9.0,9.4],timeconstants_s2[1:,i]['amp'],'bs-')
#    tosave =array([[7.0,7.8,8.4,9.0,9.4]])
#    print tosave.shape
#    for i in info:
#        print i.shape
#        tosave = append(tosave, transpose(i[1:]), axis=0)#,flatten(timeconstants_s1), flatten(amplitudes_s1), *timeconstants_s2, *amplitudes_s2]).reshape((6,-1))
#    savetxt('160209/ratesandtimeconstants.csv', tosave, delimiter = ',', header = 'pH, t_sample1')
    return None
    
def Feb8Alec():
    """Raman Spectra of Perovskite quantum dots with and without thiophenol"""
    ax = figure().add_subplot(111)
   
    a = RamanSpectrum('160208/160208_02.txt')
    b = RamanSpectrum('160208/160208_03.txt')
    c = RamanSpectrum('160208/160208_04.txt')
    d = RamanSpectrum('160208/160208_05.txt')
    CdTP = RamanSpectrum('150507/150507_01.txt')
    CdTP.normalize()
    (MeOTPRef/max(MeOTPRef)).plot()
    (CdTP/0.4+1.2).plot()
    
    offset = 2.4
    for i in [b,c]:
        i.autobaseline((875,905),order=1,join='start')
        i.autobaseline((905,1140),order=2,join='start')
        i.autobaseline((1140,1172),order=1,join='start')
        i.autobaseline((1172,1396),order=2,join='start')
        i.autobaseline((1396,1422),order=1,join='start')
        i.autobaseline((1422,1865),order=2,join='start')
        i.autobaseline((1865,1893),order=1,join='start')
        i.autobaseline((1893,2087),order=2,join='start')
        i.autobaseline((2087,2109),order=1,join='start')
        i.autobaseline((2109,2296),order=2,join='start')
        i.autobaseline((2296,2318),order=1,join='start')
        i.autobaseline((2318,2494),order=2,join='start')
        i.autobaseline((2494,2516),order=1,join='start')
       
        i.autobaseline((604,797,961,1253,1335,1494,1614,1685,1865),order=5,specialoption='points',join='start')
        
        i.autobaseline((1865,2516),order=2,join='start')
        i.autobaseline((2516,3000,3161,3605),order=2,specialoption='points',join='start')
       # i.smooth()
        i=RamanSpectrum(i.truncate(after=3590))
        i.normalize()
        i-=min(i)
        i.values[:]+=offset
        offset=max(i[600:3400])+0.2
        
        i.plot()
    
    legend(['MeOTP','Cd(thiophenol$_2$)', 'perovskite', 'perovskite+thiophenol'])
    
    CdTP.to_csv('/home/chris/Desktop/Alec/Cd-thiophenol-complex.csv',header = 'freq (cm-1), intensity')
    (MeOTPRef/max(MeOTPRef)).to_csv('/home/chris/Desktop/Alec/free_methylthiophenol.csv',header = 'freq (cm-1), intensity')
    b.to_csv('/home/chris/Desktop/Alec/perovskite_no_thiophenol.csv',  header = 'freq (cm-1), intensity')
    c.to_csv('/home/chris/Desktop/Alec/perovskite_with_thiophenol.csv', header = 'freq (cm-1), intensity')
    ylabel('Intensity')
    xlabel('Frequency (cm-1)')
    xlim(100,3600)
    ylim(0,6.5)
    savefig('/home/chris/Desktop/Alec/perovskite_figure.png', dpi=300)
    
    
    return 0
    



def Feb9UVVis(): 
    """UVVis of CdS quantum dots used for TCSPC on FEb9"""
    uvvisfile = '160209/160209_TCSPCdots.csv'
    a = loadtxt(uvvisfile, delimiter = ',', unpack =True, skiprows =1)
    oleatepeak = findpeak(a[0],a[11],(405,412))
    print 'oleate peak:', oleatepeak
    peaks = array([])
    for spec in a[1:-2]:
       # plot(a[0],spec)
        peak = findpeak(a[0],spec, (400,415))
        peaks = append(peaks, peak[0])
    peaks = peaks.reshape((-1,2))
    pHs =[7.0,7.8,8.4,9.0,9.4]
    plot(pHs,1240000/peaks[:,0]-1240000/oleatepeak[0],'s-')
    plot(pHs,1240000/peaks[:,1]-1240000/oleatepeak[0],'s-')
    return 0
    
def Feb9QY():
    
    """calculate quantum yields of PPA capped QDs measured on Feb9 with TCSPC.  There is a correction in here, because
    I did not have a fluorescence standard on Feb9, so I will use the fluorescence of the QDs in hexanes to compare between 
    Feb10 and Feb9."""
    print '--------------calculating correction factor between Feb9 and Feb 10------------'
    uvvisfile = '160209/160209fluor/DataForQY/Oleateandanthraceneuvvis.csv'
    
    
    ax1=subplot(221)
    ax2=subplot(222)
    ax3=subplot(223)
    ax4=subplot(224)
    
    fig2 = figure()
    ax5 = fig2.add_subplot(111)
    
    ####### I used the quantum dot fluorescence spectrum from yesterday to relate today's anthracene standard to yesterday's
    ###### spectra.  No this is not ideal, but it will introduce only a constant multiplicative factor into all the 
    ##### QYs.  It will not affect the ratio of k_r's at all.  
    QYdots_Feb10 = indivQY(uvvisfile, 'b','c','160209/160209fluor/DataForQY/oleatedots_standardize.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (382,467), excitationwavelength=350,nliq=1.375)#,UVVisplot = ax1,fluorplot = ax1)
    
    print 'QY of oleate capped QDs', QYdots_Feb10 
    a = loadtxt('160209/160209fluor/DataForQY/Oleateandanthraceneuvvis.csv', unpack = True, delimiter = ',', skiprows=1, usecols=(0,1))  
    
    print 'absorbance of oleate capped dots at 350:', a[1,argmin(abs(a[0]-350))],'; absorbance at 375:', a[1,argmin(abs(a[0]-375))] 
    oleateabsfactor = (1-10**-a[1,argmin(abs(a[0]-350))])/(1-10**-a[1,argmin(abs(a[0]-375))])
    print 'oleate absfactor', oleateabsfactor   ## correction for absorbance at 375 and 350.  Fluorescence of CdS on Feb9, excited at 375.  That on Feb10 excited at 350nm)
   
   
   
    QYdots_Feb9 = oleateabsfactor*indivQY(uvvisfile, 'b','c','160209/160209fluor/CdSOleatehex.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',subtractfluorfile =None, fluorescencerange = (382,467), nliq=1.375, excitationwavelength=350,fluorplot = ax2,_plot_standard=True)
    
    multfactor = QYdots_Feb10/QYdots_Feb9
    print 
    print 'mulitply all valued of fluoresence yield by', multfactor 
    
   
#    
#    """check the QY correction by commparing QY for CdS2 at pH8 from yesterday and today.  There must be a correction for the difference in excitation wavelength"""
#    a = loadtxt('160209/160209_TCSPCdots.csv', unpack = True, delimiter = ',', skiprows=1, usecols=(0,4))  
#    
#    print 'absorbance at 350 on Feb 9:', a[1,argmin(abs(a[0]-350))],'; absorbance at 375 on Feb10:', a[1,argmin(abs(a[0]-375))] 
#    absfactor = a[1,argmin(abs(a[0]-350))]/a[1,argmin(abs(a[0]-375))]  
#    print "Feb 9 spectrum should be multiplied by factor of", absfactor
#    
#    early = (absfactor*multfactor*indivQY('160209/160209_TCSPCdots.csv', 'e','m','160209/160209fluor/early_CdS2pH8.dat', 
#                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467),subtractfluorfile =None, excitationwavelength=350,UVVisplot = None,fluorplot = ax5))
#    '160209/160209fluor/BeforeTCSPC/water1.dat'
#    late =multfactor*indivQY('160209/160209_TCSPCdots.csv', 'd','m','160209/160209fluor/CdS1pH8.dat', 
#                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = None,fluorplot = ax5)
#    print early, late
#    return 0
    ###########################3
    sample1 = list()
    sample2 = list()
    uvvisfile = '160209/160209_TCSPCdots.csv'
    sample1.append(multfactor*indivQY(uvvisfile, 'b','m','160209/160209fluor/CdS1pH6.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax1,fluorplot = ax2))
    
    sample1.append(multfactor*indivQY(uvvisfile, 'd','m','160209/160209fluor/CdS1pH8.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(multfactor*indivQY(uvvisfile, 'f','m','160209/160209fluor/CdS1pH9.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(multfactor*indivQY(uvvisfile, 'h','m','160209/160209fluor/CdS1pH10.dat',
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(multfactor*indivQY(uvvisfile, 'j','m','160209/160209fluor/CdS1pH11.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax1,fluorplot = ax2))
    
    ax1.legend(['6','8','9','10','11'])
    
   

    sample2.append(multfactor*indivQY(uvvisfile, 'c','m','160209/160209fluor/CdS2pH6.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax3,fluorplot = ax4))
    
    sample2.append(multfactor*indivQY(uvvisfile, 'e','m','160209/160209fluor/CdS2pH8.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax3,fluorplot = ax4))
    sample2.append(multfactor*indivQY(uvvisfile, 'g','m','160209/160209fluor/CdS2pH9.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax3,fluorplot = ax4))
    sample2.append(multfactor*indivQY(uvvisfile, 'i','m','160209/160209fluor/CdS2pH10.dat',
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax3,fluorplot = ax4))
    sample2.append(multfactor*indivQY(uvvisfile, 'k','m','160209/160209fluor/CdS2pH11.dat', 
                    '160209/160209fluor/DataForQY/anthracene.dat',fluorescencerange = (404,467), excitationwavelength=350,UVVisplot = ax3,fluorplot = ax4))
    
    ax3.legend(['6','8','9','10','11'])
    
    return sample1,sample2

    
def Feb9_calculate_kr():
    """calculated radiative rate of recombination (kr) for Feb9 data.  TCSPC and Fluorescence"""
    plt.close('all')    
    pHs =[7.0,7.8,8.4,9.0,9.4]
    QY1,QY2 = Feb9QY()
#    Feb9TCSPC(3)
#    
#    tau_ave_s1 = sum(timeconstants_s1['tau']**2*timeconstants_s1['amp'],axis = 1)/sum(timeconstants_s1['tau']*timeconstants_s1['amp'],axis = 1)
#    tau_ave_s2 = sum(timeconstants_s2['tau']**2*timeconstants_s2['amp'],axis = 1)/sum(timeconstants_s2['tau']*timeconstants_s2['amp'],axis = 1)
#    print tau_ave_s1, tau_ave_s2, QY1, QY2
#    
#    
#    kr_s1_3exp = array(QY1)/tau_ave_s1[1:]
#    kr_s2_3exp = array(QY2)/tau_ave_s2[1:]
    
    (timeconstants_s1, timeconstants_s2) = Feb9TCSPC(2)
    tau_ave_s1 = sum(timeconstants_s1['tau']**2*timeconstants_s1['amp'],axis = 1)/sum(timeconstants_s1['tau']*timeconstants_s1['amp'],axis = 1)
    tau_ave_s2 = sum(timeconstants_s2['tau']**2*timeconstants_s2['amp'],axis = 1)/sum(timeconstants_s2['tau']*timeconstants_s2['amp'],axis = 1)
    print tau_ave_s1, tau_ave_s2, QY1, QY2
 
    
    kr_s1_2exp = array(QY1)/tau_ave_s1[1:]
    kr_s2_2exp = array(QY2)/tau_ave_s2[1:]
    
#    Feb9TCSPC(1)
#    
#    tau_ave_s1 = sum(timeconstants_s1['tau']**2*timeconstants_s1['amp'],axis = 1)/sum(timeconstants_s1['tau']*timeconstants_s1['amp'],axis = 1)
#    tau_ave_s2 = sum(timeconstants_s2['tau']**2*timeconstants_s2['amp'],axis = 1)/sum(timeconstants_s2['tau']*timeconstants_s2['amp'],axis = 1)
#    print tau_ave_s1, tau_ave_s2, QY1, QY2
#    
#    
#    kr_s1_1exp = array(QY1)/tau_ave_s1[1:]
#    kr_s2_1exp = array(QY2)/tau_ave_s2[1:]
    knrs = 1/timeconstants_s1['tau'][1:]-transpose([QY1/tau_ave_s1[1:]]*2)
    knrs2 = 1/timeconstants_s2['tau'][1:]-transpose([QY2/tau_ave_s2[1:]]*2)
    close('all')
    figure()
    suptitle('non radiative rates')
    plot(pHs, knrs[:,0])
    plot(pHs, knrs[:,1])
    plot(pHs, knrs2[:,0])
    plot(pHs, knrs2[:,1])
    legend(['sample1pop1', 'sample1pop2', 'sample2pop1', 'sample2pop2'])
    
    figure()
    

#    plot(pHs, kr_s1_3exp,'r')
#    plot(pHs, kr_s2_3exp,'r')
    
    plot(pHs, kr_s1_2exp,'b')
    plot(pHs, kr_s2_2exp,'b')
    
#    plot(pHs, kr_s1_1exp,'k')
#    plot(pHs, kr_s2_1exp,'k')
    xlabel('pH')
    ylabel('k$_r$ (ns$^{-1}$)')
    
    QYofOleateCappedQDs = 0.0697
    kr_oleate = QYofOleateCappedQDs/tau_ave_s1[0]
    print '---------final results----------'
    print 'k$_r$ of oleate capped dots', kr_oleate
#    print 'kr_s1_exp3:', kr_s1_3exp,'\n'
#    print 'kr_s2_exp3:', kr_s2_3exp,'\n'
    print 'kr_s1_exp2:', kr_s1_2exp,'\n'
    print 'kr_s2_exp2:', kr_s2_2exp,'\n'

    
#    
#    figure()
#    plot(pHs,tau_ave_s1[1]/tau_ave_s1[1:],'--')
#    plot(pHs, tau_ave_s2[1]/tau_ave_s2[1:],'--')
#    plot(pHs, QY1/QY1[0])
#    plot(pHs, QY2/QY2[0])
#    plot(pHs, kr_s1_2exp/kr_s1_2exp[0],'b')
#    plot(pHs, kr_s2_2exp/kr_s2_2exp[0],'b')
    
    return (np.append(pHs, pHs), np.append(kr_s1_2exp,kr_s2_2exp))
    
def Feb16TCSPC(numberofexponentials, axis):
    """TCSPC results from PPA capped CdS dots on Feb9 at pH 6 8 and 11"""
   # os.chdir('/media/cybertron_box/Chris Thompson/TCSPC Data')

    files = ['160216/2-16-16-375exc--irf-4mhz-250psbin-750em-scattercell',
             '160216/2-16-16-375exc-cdsoleate-hexanes-chris-522p602t-250psbin-1496rate-1mhz',
             '160216/2-16-16-375exc-cdsph11-chris-273p157t-250psbin-297rate-8mhz',
             '160216/2-16-16-375exc-cdsph7-s2-chris-505p96t-250psbin-281rate-4mhz',
             '160216/2-16-16-375exc-cdsph11-s2-chris-2312p9t-100psbin-300rate-4mhz',
             '160216/2-16-16-375exc-cdsph5-s2-chris-883p1t-250psbin-311rate-4mhz',
             '160216/2-16-16-375exc-cdsph10-s2-chris-270p1t-250psbin-309rate-4mhz',
             '160216/2-16-16-375exc-cdsph12-s2-chris-3423p8t-250psbin-308rate-4mhz',
            '160216/2-16-16-375exc-cdsph5-s1-chris-2263p2t-250psbin-309rate-4mhz',
            '160216/2-16-16-375exc-cdsph7-s1-chris-1638p3t-250psbin-305rate-4mhz',
            '160216/2-16-16-375exc-cdsph10-s1-chris-954p3t-250psbin-305rate-4mhz',
            '160216/2-16-16-375exc-cdsph8p5-s2-chris-192p4t-250psbin-329rate-4mhz',
            '160216/2-16-16-375exc-cdsph8p5-s1-chris-434p8t-250psbin-310rate-4mhz',
            '160216/2-16-16-375exc-cdsph12-s2-chris-3681p25t-250psbin-304rate-4mhz'            ,
            ]
    IRF_1 = loadtxt(files[0],unpack = True)
    oleate_1 = loadtxt(files[1],unpack = True,)

    CdS1pH11 = loadtxt(files[2],unpack = True)
    CdS2pH7 = loadtxt(files[3],unpack = True)
    CdS2pH11 = loadtxt(files[4],unpack = True)
    CdS2pH5 = loadtxt(files[5],unpack = True)
    CdS2pH10 = loadtxt(files[6],unpack = True)
    CdS2pH12 = loadtxt(files[7],unpack = True)
    CdS1pH5 = loadtxt(files[8],unpack = True)
    CdS1pH7 = loadtxt(files[9],unpack = True)
    CdS1pH10 = loadtxt(files[10],unpack = True)
    CdS2pH9= loadtxt(files[11],unpack = True)
    CdS1pH9= loadtxt(files[12],unpack = True)
    CdS1pH12 =loadtxt(files[13],unpack = True)
    
   
    
    
    for i  in [IRF_1,oleate_1,CdS1pH11,CdS1pH5,CdS1pH12, CdS1pH10, CdS1pH7,CdS2pH5,CdS2pH7,CdS2pH10,CdS2pH11,CdS2pH12,CdS2pH9, CdS1pH9]:#, CdS2pH6, CdS2pH8,CdS2pH9, CdS2pH10, CdS2pH11]:
       # i[1] = SGsmooth(i[0],i[1])
        zero = numpy.mean(i[1,10:argmax(i[1])-10],axis = 0)
     
        i[1]-=zero
        i[1]/=max(i[1])
        
 


    """FIT THE IRF TO A GAUSSIAN"""
    def gaussian(x,A,w,x0):return (A/sqrt(2*pi)/w)*exp(-(x-x0)**2/(2*w**2))
        
    fitIRF = scipy.optimize.curve_fit(gaussian, IRF_1[0], IRF_1[1],[1,1,45])[0]
    #plot(IRF_1[0], gaussian(IRF_1[0],*fitIRF))
    fixt0=fitIRF[2]
    
    w0 = fitIRF[1]

    
    print 'IRFfits', w0, fixt0
    figure()
    plot(IRF_1[0], IRF_1[1])
    plot(IRF_1[0], gaussian(IRF_1[0],*fitIRF))
 

    
    ## FIT THE DATA TO THE IRF CONVOLUTED WITH EXPONENTIALS
    def fitfunction(x,*args):
        """t0, A1, tau1, A2, tau2..."""
        y = numpy.zeros(x.shape)
        t0=args[0]
        for z in range((len(args)-1)/2):
            
            A1 = args[z*2+1]
            tau1 = args[z*2+2]
           
            y+=(A1/2)*exp(w0**2/(2*tau1**2) - (x-t0)/tau1)* (1 - erf((w0**2-tau1*(x-t0))/(w0*tau1*sqrt(2))))   # / (exp(w0**2/(2*tau1**2))* (1 - erf((w0)/(tau1*sqrt(2)))))
        return y
    def fitfunction2(x,t0, A1, tau1,A2):
        """t0, A1, tau1, A2, tau2..."""

           
        y=(A1/2)*exp(w0**2/(2*tau1**2) - (x-t0)/tau1)* (1 - erf((w0**2-tau1*(x-t0))/(w0*tau1*sqrt(2)))) +(A2/2)*exp(w0**2/(2*0.189**2) - (x-t0)/0.189)* (1 - erf((w0**2-0.189*(x-t0))/(w0*0.189*sqrt(2))))   # / (exp(w0**2/(2*tau1**2))* (1 - erf((w0)/(tau1*sqrt(2)))))
        return y
    figure()
    suptitle('sample1')
    
 
  
    errors = array([])
    timeconstants_s1 = ndarray((0,3),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        
    print '-----performing fits -----'
    for r  in range(1,8):
        
        subplot(710+r)
        if numberofexponentials==3:
            if r==0:
                guess = [44,0.05,12,0.2,6,0.04,0.5]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6,0.2,1,0.04,0.5]#,0.04, 1]
            else:
                guess = [44.3,0.05,40, 0.01,25,    0.06,   5]
        elif numberofexponentials ==2:
            if r==0:
                guess = [44,0.05,12,0.04,0.5]#,0.04, 1]
            else:
                guess = [44,0.5,6,0.5,0.5]#,0.04, 1]
            
        elif numberofexponentials ==1:
            
            guess = [44,1,1]#,0.04, 1]
            
        data = list([oleate_1, CdS1pH5,CdS1pH7,CdS1pH9,CdS1pH10,CdS1pH11,CdS1pH12 ])[r-1]
        if type(data) is int:
            continue

        #print numpy.mean(numpy.diff(i[0]))
        plot(data[0],data[1],'.')  #### plot the data
        
       # plot(i[0],fitfunction(i[0],*guess))   ### plot the guess
        try:
            fit = scipy.optimize.curve_fit(fitfunction, data[0], data[1],guess,maxfev=3000)[0]
            
        except:
            fit = array(guess)
        #v=fit[1] + fit[3] +fit[5]
        xlim(40,60)
        ylim(-0.05,1.05)
       
        print list(("%.3f" % i for i in fit))
        plot(data[0],fitfunction(data[0],*fit),linewidth=3) ### plot the total fit
        errorstart = argmin(abs(42-data[0]))
        errors = np.append(errors,sqrt(sum(((data[1][errorstart:]-fitfunction(data[0][errorstart:],*fit))/data[1][errorstart:])**2)))
        
        # ratearray = array([(fit[1], fit[2]), (fit[3], fit[4]), (fit[5], fit[6])],     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray = array(list((fit[i+1], fit[i+2]) for i in arange(numberofexponentials)*2),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray=np.sort(ratearray,order = 'tau')
        for i in ratearray:
            plot(data[0],fitfunction(data[0],fit[0],i['amp'],i['tau']),linewidth=1)
    
        

        timeconstants_s1= np.append(timeconstants_s1,ratearray)
        
        annotate(['oleate','PPA pH5','PPA pH7','pH8.5', 'pH10','pH11','pH12'][r-1],(0.8,0.8),xycoords='axes fraction')
    timeconstants_s1=timeconstants_s1.reshape((-1,numberofexponentials))
    
    xlabel('time (ns)')
    

    timeconstants_s2 =ndarray((0,3),     dtype=[('amp', '<f8'), ('tau', '<f8')])
    
    figure()
    suptitle('sample2')
    for r  in range(1,8):
        
        subplot(710+r)
        if numberofexponentials==3:
            if r==0:
                guess = [44,0.05,12,0.2,6,0.04,0.5]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6,0.2,1,0.04,0.5]#,0.04, 1]
            else:
                guess = [44.3,0.05,40, 0.01,25,    0.06,   5]
        elif numberofexponentials ==2:
            if r==0:
                guess = [44,0.05,12,0.04,0.5]#,0.04, 1]

            else:
                guess = [44.3,0.05,40,    0.06,   5]
        elif numberofexponentials ==1:
            if r==0:
                guess = [44,0.05,12]#,0.04, 1]
            elif r==6:
                guess = [44,0.05,6]#,0.04, 1]
            else:
                guess = [44.3,0.05,40]
        data = list([oleate_1, CdS2pH5, CdS2pH7,CdS2pH9,CdS2pH10, CdS2pH11, CdS2pH12])[r-1]
        #print numpy.mean(numpy.diff(i[0]))
        plot(data[0],data[1],'.')  #### plot the data
       # plot(i[0],fitfunction(i[0],*guess))   ### plot the guess
        try:
            fit = scipy.optimize.curve_fit(fitfunction, data[0], data[1],guess,maxfev=3000)[0]
        except:
            fit = array(guess)
        #v=fit[1] + fit[3] +fit[5]
        xlim(40,50)
        ylim(-0.05,1.05)
       
        print list(("%.3f" % i for i in fit))
        plot(data[0],fitfunction(data[0],*fit),linewidth=3) ### plot the total fit
        errorstart = argmin(abs(42-data[0]))
        errors = np.append(errors,sqrt(sum(((data[1][errorstart:]-fitfunction(data[0][errorstart:],*fit))/data[1][errorstart:])**2)))
        
        
        ratearray = array(list((fit[i+1], fit[i+2]) for i in arange(numberofexponentials)*2),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray=np.sort(ratearray,order = 'tau')
        for i in ratearray:
            plot(data[0],fitfunction(data[0],fit[0],i['amp'],i['tau']),linewidth=1)
       
        

        timeconstants_s2= np.append(timeconstants_s2,ratearray)
        annotate(['oleate','PPA pH6','PPA pH8','pH9','pH10','PPA pH11','a'][r-1],(0.8,0.8),xycoords='axes fraction')
    timeconstants_s2=timeconstants_s2.reshape((-1,numberofexponentials))
    
    print 'average error', np.mean(errors)
    xlabel('time (ns)')
    
    return (timeconstants_s1, timeconstants_s2)
    
def compareFeb16andFeb9():
    a = Feb16TCSPC(2,gca())
    b = Feb9TCSPC(2)
    c = array(list(sum(i['tau']*i['amp'])/sum(i['tau']**2*i['amp']) for i in a[0]))
    d = array(list(sum(i['tau']*i['amp'])/sum(i['tau']**2*i['amp']) for i in a[1]))
    e = array(list(sum(i['tau']*i['amp'])/sum(i['tau']**2*i['amp']) for i in b[0]))
    f = array(list(sum(i['tau']*i['amp'])/sum(i['tau']**2*i['amp']) for i in b[1]))
 

   
    kr16_1 = array([  5.68940824e-06,   6.96327361e-06,   3.20817230e-05,
         2.61565900e-05,   1.23618275e-05,   6.01393212e-06])
    kr16_2 = array([  8.60710635e-06,   1.18440562e-05,   4.57037182e-05,
         9.00143630e-05,   1.88734351e-05,   8.28018988e-07])
    
    kr_s1_exp2 = array([  2.38746485e-05,   3.87144561e-05,   6.44469412e-05,
                        9.53312634e-05,   2.19309118e-04] )

    kr_s2_exp2 =array([  2.01436407e-05 ,  7.42887502e-05 ,  6.94277829e-05 ,  1.59649983e-04,
                 2.94667399e-04])
   
    close('all')
    figure()
    plot([5.6, 7.5,8.4,10.0,11.1,11.8],c[1:],'ks')
    plot([5.6, 7.5,8.4,10.2,11.1,12.0],d[1:],'ks')
    plot([7.0,7.8,8.4,9.0,9.4],e[1:],'ks')
    plot([7.0,7.8,8.4,9.0,9.4],f[1:],'ks')
    
    figure()
    pHsaggregated = array([5.6, 7.5,8.4,10.0,11.1,11.8,5.6, 7.5,8.4,10.0,11.1,11.8,7.0,7.8,8.4,9.0,9.4,7.0,7.8,8.4,9.0,9.4])
    QYsaggregated = append(kr16_1/kr16_1[2],kr16_2/kr16_2[2])
    QYsaggregated= append(QYsaggregated,kr_s1_exp2/kr_s1_exp2[2])
    QYsaggregated = append(QYsaggregated,kr_s2_exp2/kr_s2_exp2[2])
    plot(pHsaggregated,QYsaggregated,'sk')
    fitted = numpy.polynomial.polynomial.polyfit(pHsaggregated,QYsaggregated,5)
    plot(linspace(5.6,12,100),numpy.polynomial.polynomial.polyval(linspace(5.6,12,100),fitted),'r')
#    plot([5.6, 7.5,8.4,10.0,11.1,11.8],kr16_1/kr16_1[2],'sk')
#    plot([5.6, 7.5,8.4,10.2,11.1,12.0],kr16_2/kr16_2[2],'sk')
#    plot([7.0,7.8,8.4,9.0,9.4],kr_s1_exp2/kr_s1_exp2[2],'sk')
#    plot([7.0,7.8,8.4,9.0,9.4],kr_s2_exp2/kr_s2_exp2[2],'sk')
    
    ylabel('k$_r$ (ns$^{-1}$)')
    xlabel('pH')
    
    QYFeb91= array([0.00010049827821621676, 0.00013066903216727769, 0.00022781724298307344, 0.00025977476507963978, 0.00023207098018508698])*100
    QYFeb92=array([8.4493517569121898e-05, 0.00015499338993774661, 0.00023270709386317347, 0.00030337966297159114, 0.00026606606002984402])*100
    QYFeb161 = array([  1.22096847e-05,
                     1.63693716e-05,   2.63285990e-05,
         1.93740716e-05,   6.82283240e-06,   6.99873779e-06])*100
    QYFeb162 = array([  2.64292285e-05,   3.57962025e-05,   7.18266873e-05,
         4.94787370e-05,   1.16480806e-05,   8.31632682e-07])*100
    figure() 
    
    plot([5.6, 7.5,8.4,10.0,11.1,11.8],QYFeb161,'-s')
    plot([5.6, 7.5,8.4,10.2,11.1,12.0],QYFeb162,'-s')
    plot([7.0,7.8,8.4,9.0,9.4],QYFeb91,'-s')
    plot([7.0,7.8,8.4,9.0,9.4],QYFeb92,'-s')
    ylabel('PL QY (%)')
    xlabel('pH')
    legend(['Day1_sample1', 'Day1_sample2', 'Day2_sample1', 'Day2_sample2'])
    return 0
def Feb16UVVis():
    uv = loadtxt('160216/CdSPPA_TCSPC_160216.csv',skiprows = 2,delimiter = ',', unpack = True,usecols = (0,1,3,12,5,7,9))
    uv2 = loadtxt('160216/CdSPPA_TCSPC_160216.csv',skiprows = 2,delimiter = ',', unpack = True,usecols = (0,2,4,13,6,8,10))
    uvoleate  = loadtxt('160216/CdSPPA_TCSPC_160216.csv',skiprows = 2,delimiter = ',', unpack = True,usecols = (0,15))
    
    oleatepeak = findpeak(uvoleate[0],uvoleate[1],(405,415))
    peaks1 = array([])
    peaks2 = array([])
    for i in uv[1:]:
        plot(uv[0],i)
        r = findpeak(uv[0],i,(405,425))
        print r
        peaks1 = append(peaks1,1240000/r[0]-1240000/oleatepeak[0])
    plot(uvoleate[0],uvoleate[1])
    for i in uv2[1:]:
         
        
        r = findpeak(uv2[0],i,(405,425))
        peaks2 = append(peaks2,1240000/r[0]-1240000/oleatepeak[0])
    print peaks1.shape
    
    
    fpeaks1 = array([])
    filenamelist = ['CdS1pH5.dat',
                'CdS1pH7.dat',
                'CdS1pH8pt5b.dat',
                'CdS1pH10.dat',
                'CdS1pH11.dat',
                'CdS1pH12.dat',]
    a = loadtxt('160216/'+'oleate1nm.dat', unpack=True,skiprows = 2, usecols = (0,3))
        
        
    r = findpeak(uv2[0],i,(420,450))
    oleatepeak = r[0]       
    print oleatepeak
                
    for name in filenamelist:
        a = loadtxt('160216/'+name, unpack=True,skiprows = 2, usecols = (0,3))
        
        r = findpeak(uv2[0],i,(430,450))
        fpeaks1 = append(fpeaks1,1240000/r[0]-1240000/430)
   
   
                
    filenamelist = ['CdS2pH5.dat',
                    'CdS2pH7.dat',
                    'CdS2pH8pt5.dat',
                    'CdS2pH10.dat',
                    'CdS2pH11.dat',
                    'CdS2pH12.dat',]
    for name in filenamelist:
        a = loadtxt('160216/'+name, unpack=True,skiprows = 2, usecols = (0,3))
        
        plot(a[0],a[1])
        r = findpeak(uv2[0],i,(430,450))
        fpeaks1 = append(fpeaks1,1240000/r[0]-1240000/430)                
    figure()
    plot([5.6, 7.5,8.4,10.0,11.1,11.8],(peaks1+peaks2)/2,'ks-')
    
    xlabel('pH',fontsize =20)
    ylabel('$\Delta\lambda$ relative to oleate (eV)',fontsize = 18)
    fpeaks1 = 1240000/array([421,422.5,427,426.5,432.4,443.3],dtype=float)-1240000/oleatepeak
    fpeaks2 = 1240000/array( [421,424,427.7,428.4,434,439.3],dtype=float)-1240000/oleatepeak
  
    plot([5.6, 7.5,8.4,10.0,11.1,11.8],(fpeaks1+fpeaks2)/2,'rs-')
    legend(['absorption','photoluminescence'])
    
    return 0
    
def Feb16QY():
    
    nE = 1.359
    nQ = 1.44
    nW = 1.333  ## refractive index water    
    nhex = 1.375
    
    fig = figure()
    ax1 = fig.add_subplot(111)
    ax1.set_ylabel('relative QY')
    ax1.set_xlabel('pH')
   
    filenamelist = ['oleate1nm.dat','oleate2nm.dat','anthracene1.dat','anthracene2.dat','anthracene3.dat','anthracene5.dat']
    uv = loadtxt('160216/CdSPPA_TCSPC_160216.csv',skiprows = 2,delimiter = ',', unpack = True,
                 usecols = (0,15,15,14,14,14,14))
    uv[1:]-=transpose([uv[1:,0]])
    x330=argmin(abs(uv[0]-330))
    absorbances = 1-10**-uv[1:,x330]
    print absorbances
    figure()
    for i in uv[1:]:plot(uv[0],i)
    legend(filenamelist)
        
    
    figure()    
    QYs = array([])
    for i in range(len(filenamelist)):
        
        a = loadtxt('160216/'+filenamelist[i], skiprows=2, delimiter = '\t', unpack = True, usecols = (2,3))
        
        a[1]/=absorbances[i]
        
        r = RamanSpectrum(pandas.Series(a[1],a[0]))
        if i>1:
            r-=r[550]
        else:
            r-=r[355]
        QYs=append(QYs,r[470]*1257)#r.calc_area((335,500)))
        if i ==5:
            anthracenefluorescencearea = r[470]*1257
        print r.calc_area((335,500))/r[470]
        r.plot(label=filenamelist[i])
    
    legend(filenamelist)
    oleatevalue = QYs[0]
    print 'oleate QY', oleatevalue*0.27/(1+0.00145*158)*nhex**2/nE**2/QYs[2]
    
   
    
    
    filenamelist = ['CdS1pH5.dat',
                    'CdS1pH7.dat',
                    'CdS1pH8pt5b.dat',
                    'CdS1pH10.dat',
                    'CdS1pH11.dat',
                    'CdS1pH12.dat',]
    uv = loadtxt('160216/CdSPPA_TCSPC_160216.csv',skiprows = 2,delimiter = ',', unpack = True,usecols = (0,1,3,12,5,7,9))

    uv[1:]-=transpose([uv[1:,0]])
    x330=argmin(abs(uv[0]-330))
    absorbances = 1-10**-uv[1:,x330]
    print absorbances
    figure()
    for i in uv[1:]:plot(uv[0],i)
    
    figure()
    QYs_1 = array([])
    for i in range(len(filenamelist)):
        
        a = loadtxt('160216/'+filenamelist[i], skiprows=2, delimiter = '\t', unpack = True, usecols = (2,3))
        a[1]/=absorbances[i]
        a[1]*=0.27/(1+0.00145*158)*nW**2/nE**2 /anthracenefluorescencearea
        r = RamanSpectrum(pandas.Series(a[1],a[0]))
        r-=r[400]
        QYs_1=append(QYs_1,r.calc_area((400,460)))
        r.plot(label=filenamelist[i])
    legend()
    figure()
    ax1.plot([5.6, 7.5,8.4,10.0,11.1,11.8], QYs_1)
    
    
    filenamelist = ['CdS2pH5.dat',
                    'CdS2pH7.dat',
                    'CdS2pH8pt5.dat',
                    'CdS2pH10.dat',
                    'CdS2pH11.dat',
                    'CdS2pH12.dat',]
    uv = loadtxt('160216/CdSPPA_TCSPC_160216.csv',skiprows = 2,delimiter = ',', unpack = True,usecols = (0,2,4,13,6,8,10))
    uv[1:]-=transpose([uv[1:,0]])
    x330=argmin(abs(uv[0]-330))
    absorbances = 1-10**-uv[1:,x330]
    print absorbances
    
    QYs_2 = array([])
    for i in range(len(filenamelist)):
        
        a = loadtxt('160216/'+filenamelist[i], skiprows=2, delimiter = '\t', unpack = True, usecols = (2,3))
        a[1]/=absorbances[i]
        a[1]*=0.27/(1+0.00145*158)*nW**2/nE**2 /anthracenefluorescencearea
        r = RamanSpectrum(pandas.Series(a[1],a[0]))
        r-=r[400]
        QYs_2 =append(QYs_2 ,r.calc_area((400,460)))
        r.plot(label=filenamelist[i])
      
    legend()
    figure()
    for i in uv[1:]:plot(uv[0],i)
    ax1.plot([5.6, 7.5,8.4,10.0,11.1,11.8], QYs_2)
        
            
    return (QYs_1,QYs_2)
    
def Feb16QY2():
    fig =figure()
    
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    sample1 = list()
    sample2 = list()
    


    
    
    uvvisfile = '160216/CdSPPA_TCSPC_160216.csv'
    
    
    sample1.append(indivQY(uvvisfile, 'b','o','160216/CdS1pH5.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(indivQY(uvvisfile, 'd','o','160216/CdS1pH7.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(indivQY(uvvisfile, 'm','o','160216/CdS1pH8pt5b.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(indivQY(uvvisfile, 'f','o','160216/CdS1pH10.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(indivQY(uvvisfile, 'h','o','160216/CdS1pH11.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample1.append(indivQY(uvvisfile, 'j','o','160216/CdS1pH12.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample2.append(indivQY(uvvisfile, 'c','o','160216/CdS2pH5.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample2.append(indivQY(uvvisfile, 'e','o','160216/CdS2pH7.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample2.append(indivQY(uvvisfile, 'n','o','160216/CdS2pH8pt5.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample2.append(indivQY(uvvisfile, 'g','o','160216/CdS2pH10.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample2.append(indivQY(uvvisfile, 'i','o','160216/CdS2pH11.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    sample2.append(indivQY(uvvisfile, 'k','o','160216/CdS2pH12.dat', 
                    '160216/anthracene5.dat',fluorescencerange = (400,460), excitationwavelength=330,UVVisplot = ax1,fluorplot = ax2))
    figure()
    suptitle('qy2')
    plot([5.6, 7.5,8.4,10.0,11.1,11.8],sample1)
    plot([5.6, 7.5,8.4,10.0,11.1,11.8],sample2)
    return sample1,sample2
def Feb16_calculate_kr():
    """calculated radiative rate of recombination (kr) for Feb9 data.  TCSPC and Fluorescence"""
    plt.close('all')    
    pHs =[5.6, 7.5,8.4,10.0,11.1,11.8]
    QY1,QY2 = Feb16QY()

    
    (timeconstants_s1, timeconstants_s2) = Feb16TCSPC(2,None)
    tau_ave_s1 = sum(timeconstants_s1['tau']**2*timeconstants_s1['amp'],axis = 1)/sum(timeconstants_s1['tau']*timeconstants_s1['amp'],axis = 1)
    tau_ave_s2 = sum(timeconstants_s2['tau']**2*timeconstants_s2['amp'],axis = 1)/sum(timeconstants_s2['tau']*timeconstants_s2['amp'],axis = 1)
    print tau_ave_s1, tau_ave_s2, QY1, QY2
      
    
    kr_s1_2exp = array(QY1)/tau_ave_s1[1:]
    kr_s2_2exp = array(QY2)/tau_ave_s2[1:]
    
    knrs = 1/timeconstants_s1['tau'][1:]-transpose([QY1/tau_ave_s1[1:]]*2)
    knrs2 = 1/timeconstants_s2['tau'][1:]-transpose([QY2/tau_ave_s2[1:]]*2)
    close('all')
    figure()
    suptitle('non radiative rates')
    plot(pHs, knrs[:,0])
    plot(pHs, knrs[:,1])
    plot(pHs, knrs2[:,0])
    plot(pHs, knrs2[:,1])
    legend(['sample1pop1', 'sample1pop2', 'sample2pop1', 'sample2pop2'])
    
    
    
    figure()
    suptitle('final results')

    
    plot(pHs, kr_s1_2exp,'b')
    plot(pHs, kr_s2_2exp,'b')

    xlabel('pH')
    ylabel('k_r (ns^-1)')
    
    QYofOleateCappedQDs = 0.0697
    kr_oleate = QYofOleateCappedQDs/tau_ave_s1[0]
    print '---------final results----------'
    print 'k$_r$ of oleate capped dots', kr_oleate

    print 'kr_s1_exp2:', kr_s1_2exp,'\n'
    print 'kr_s2_exp2:', kr_s2_2exp,'\n'
    

    return (kr_s1_2exp,kr_s2_2exp)
    
def Feb17CdSeTCSPC(numberofexponentials):
    """TCSPC results from PPA capped CdS dots on Feb9 at pH 6 8 and 11"""
   # os.chdir('/media/cybertron_box/Chris Thompson/TCSPC Data')

    files = ['160217/2-17-16-450exc-1mhz-cdseoleate-550em-1279p8t-50psbin-392rate',
              '160217/2-17-16-450exc-2mhz-irf-scattercell-250psbin', 
              '160217/2-17-16-450exc-2mhz-cdsppa-s1-ph11-550em-77p9t-250psbin-382rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s1-ph5-550em-390p0t-250psbin-313rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s1-ph7-550em-150p4t-250psbin-304rate'
              , '160217/2-17-16-450exc-2mhz-cdsppa-s1-ph8-550em-106p6t-250psbin-343rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s1-ph9-550em-94p9t-250psbin-370rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s2-ph11-550em-113p2t-250psbin-390rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s2-ph5-550em-206p0t-250psbin-295rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s2-ph7-550em-115p5t-250psbin-337rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s2-ph8-550em-132p0t-250psbin-311rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s2-ph9-550em-91p8t-250psbin-340rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s3-ph10-550em-85p9t-250psbin-410rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s3-ph5-550em-231p5t-250psbin-296rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s3-ph7-550em-143p7t-250psbin-340rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s3-ph8-550em-144p4t-250psbin-301rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s3-ph9-550em-82p9t-250psbin-377rate', 
              '160217/2-17-16-450exc-2mhz-cdsppa-s4-ph11-550em-99p9t-250psbin-320rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s4-ph5-550em-275p0t-250psbin-274rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s4-ph7-550em-171p9t-250psbin-319rate',
              '160217/2-17-16-450exc-2mhz-cdsppa-s4-ph8-550em-182p8t-250psbin-283rate', 
              '160217/2-17-16-450exc-2mhz-cdsppa-s4-ph9-550em-122p3t-250psbin-324rate', ]
              
  
    
    ( oleate_1 ,IRF_1,
    CdS1pH5,
    CdS1pH11,
    CdS1pH7 ,
    CdS1pH8,
    CdS1pH9,
    CdS2pH5,
    CdS2pH11,
    CdS2pH7,
    CdS2pH8,
    CdS2pH9,
    CdS3pH5,
    CdS3pH11,
    CdS3pH7,
    CdS3pH8,
    CdS3pH9,
    CdS4pH5,
    CdS4pH11 ,
    CdS4pH7,
    CdS4pH8 ,
    CdS4pH9) = tuple((loadtxt(i,unpack = True) for i in files))
    
    
    
    for i  in [IRF_1, oleate_1 ,    CdS1pH5,    CdS1pH11,    CdS1pH7 ,    CdS1pH8,    CdS1pH9,    CdS2pH5,
    CdS2pH11,    CdS2pH7,    CdS2pH8,    CdS2pH9,    CdS3pH5,    CdS3pH11,    CdS3pH7,    CdS3pH8,    CdS3pH9,
    CdS4pH5,    CdS4pH11 ,    CdS4pH7,    CdS4pH8 ,    CdS4pH9]:       # i[1] = SGsmooth(i[0],i[1])
        zero = numpy.mean(i[1,10:argmax(i[1])-10],axis = 0)
     
        i[1]-=zero
        i[1]/=max(i[1])
        
 


    """FIT THE IRF TO A GAUSSIAN"""
    def gaussian(x,A,w,x0):return (A/sqrt(2*pi)/w)*exp(-(x-x0)**2/(2*w**2))
        
    fitIRF = scipy.optimize.curve_fit(gaussian, IRF_1[0], IRF_1[1],[1,1,45])[0]
    #plot(IRF_1[0], gaussian(IRF_1[0],*fitIRF))
    fixt0=fitIRF[2]
    
    w0 = fitIRF[1]

    
    print 'IRFfits', w0, fixt0
    figure()
    plot(IRF_1[0], IRF_1[1])
    plot(IRF_1[0], gaussian(IRF_1[0],*fitIRF))
 

    
    ## FIT THE DATA TO THE IRF CONVOLUTED WITH EXPONENTIALS
    def fitfunction(x,*args):
        """t0, A1, tau1, A2, tau2..."""
        y = numpy.zeros(x.shape)
        t0=args[0]
        for z in range((len(args)-1)/2):
            
            A1 = args[z*2+1]
            tau1 = args[z*2+2]
           
            y+=(A1/2)*exp(w0**2/(2*tau1**2) - (x-t0)/tau1)* (1 - erf((w0**2-tau1*(x-t0))/(w0*tau1*sqrt(2))))   # / (exp(w0**2/(2*tau1**2))* (1 - erf((w0)/(tau1*sqrt(2)))))
        return y

    
 
  
    errors = array([])
    timeconstants_s1 = ndarray((0,3),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        
    print '-----performing fits -----'
    for r  in range(1,25):
        if (r-1)%6==0:
            figure()
            suptitle('sample'+str((r+7)/6))
        
        subplot(610+(r-1)%6+1)
        if numberofexponentials==3:
            if r==0:
                guess = [fixt0,0.05,12,0.04,0.5,0.04,0.5]#,0.04, 1]
            elif r==6:
                guess = [fixt0,0.05,12,0.04,0.5,0.04,0.5]#,0.04, 1]
            else:
                guess = [fixt0,0.05,12,0.04,0.5,  0.06,   5]
        elif numberofexponentials ==2:
            if r==0:
                guess = [fixt0,0.05,12,0.04,0.5]#,0.04, 1]
            else:
                guess = [fixt0,0.5,10,0.5,1]#,0.04, 1]
            
        elif numberofexponentials ==1:
            
            guess = [fixt0,1,1]#,0.04, 1]
            
        data = list([oleate_1 ,    CdS1pH5,        CdS1pH7 ,    CdS1pH8,    CdS1pH9, CdS1pH11,
                     oleate_1 , CdS2pH5,     CdS2pH7,    CdS2pH8,    CdS2pH9,   CdS2pH11,  
                     oleate_1 , CdS3pH5,       CdS3pH7,    CdS3pH8,    CdS3pH9,CdS3pH11, 
                     oleate_1 , CdS4pH5,    CdS4pH7,    CdS4pH8 ,    CdS4pH9, CdS4pH11   ])[r-1]
        if type(data) is int:
            continue

        #print numpy.mean(numpy.diff(i[0]))
        plot(data[0],data[1],'.')  #### plot the data
        
       # plot(i[0],fitfunction(i[0],*guess))   ### plot the guess
        try:
            fit = scipy.optimize.curve_fit(fitfunction, data[0], data[1],guess,maxfev=3000)[0]
            
        except:
            fit = array(guess)
        #v=fit[1] + fit[3] +fit[5]
        xlim(40,60)
        ylim(-0.05,1.05)
       
        print list(("%.3f" % i for i in fit))
        plot(data[0],fitfunction(data[0],*fit),linewidth=3) ### plot the total fit
        errorstart = argmin(abs(42-data[0]))
        errors = np.append(errors,sqrt(sum(((data[1][errorstart:]-fitfunction(data[0][errorstart:],*fit))/data[1][errorstart:])**2)))
        
        # ratearray = array([(fit[1], fit[2]), (fit[3], fit[4]), (fit[5], fit[6])],     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray = array(list((fit[i+1], fit[i+2]) for i in arange(numberofexponentials)*2),     dtype=[('amp', '<f8'), ('tau', '<f8')])
        ratearray=np.sort(ratearray,order = 'tau')
        for i in ratearray:
            plot(data[0],fitfunction(data[0],fit[0],i['amp'],i['tau']),linewidth=1)
    
        

        timeconstants_s1= np.append(timeconstants_s1,ratearray)
        
        #annotate(['oleate','PPA pH5','PPA pH7','pH8.5', 'pH10','pH11','pH12'][r-1],(0.8,0.8),xycoords='axes fraction')
    timeconstants_s1=timeconstants_s1.reshape((-1,numberofexponentials))
    
    xlabel('time (ns)')
    

    figure()
    tau_ave = array(sum(timeconstants_s1['tau']**2*timeconstants_s1['amp'],axis=1)/sum(timeconstants_s1['tau']*timeconstants_s1['amp'],axis=1))
    tau_ave = tau_ave.reshape((4,6))
    for i in range(4):plot(tau_ave[i][1:])
#    xlabel('time (ns)')
    
    return (timeconstants_s1)
    
def Feb17UVVisfigure():
    uvvis = loadtxt('160217/160217TCSPC.csv', delimiter = ',', unpack = True,skiprows = 1,usecols=(0,1,9,13,17))
    
    for i in uvvis[1:]:
        i-=i[0]
        i=SGsmooth(uvvis[0],i)
        p = findpeak(uvvis[0],i,(520,530))
        i/=p[1]
        plot(uvvis[0],i)
    return 0
def Feb17CdSeQYs():
    fluorfilelist = ['160217/Raman160217/160217_26.txt',
                    '160217/Raman160217/160217_05.txt',
                     '160217/Raman160217/160217_04.txt',
                     '160217/Raman160217/160217_06.txt',
                     '160217/Raman160217/160217_08.txt',
                     '160217/Raman160217/160217_09.txt',
                     '160217/Raman160217/160217_10.txt',
                     '160217/Raman160217/160217_11.txt',
                     '160217/Raman160217/160217_12.txt',
                     '160217/Raman160217/160217_13.txt',
                     '160217/Raman160217/160217_14.txt',
                     '160217/Raman160217/160217_15.txt',
                     '160217/Raman160217/160217_16.txt',
                     '160217/Raman160217/160217_17.txt',
                     '160217/Raman160217/160217_18.txt',
                     '160217/Raman160217/160217_19.txt',
                     '160217/Raman160217/160217_20.txt',
                     '160217/Raman160217/160217_21.txt',
                     '160217/Raman160217/160217_22.txt',
                     '160217/Raman160217/160217_23.txt',
                     '160217/Raman160217/160217_24.txt']
                     
    filenames = array(['1pH5', '2pH5', '3pH5', '4pH5',
                       '1pH7', '2pH7', '3pH7', '4pH7',
                       '1pH8', '2pH8', '3pH8', '4pH8',
                       '1pH9', '2pH9', '3pH5', '4pH9',
                       '1pH11', '2pH11', '3pH11', '4pH11',])
                     
    uvvis = loadtxt('160217/160217TCSPC.csv', delimiter = ',', unpack = True,skiprows = 1)
    uvvis[1:]-=transpose([uvvis[1:,0]])
    forfigure = uvvis[array([0,1,4])]
    forfigure
    x473nm = np.where(uvvis[0]==473)[0][0]
    print x473nm
    
    for i in uvvis[1:]: plot(uvvis[0],i)
    figure()
    absorbances_dots = uvvis[1:-1,x473nm]
    plot(absorbances_dots)
    figure()
    absorbance_anth = uvvis[-1,x473nm]
    
    fluorescencedots = array([])
    
    for i in fluorfilelist:
        s = RamanSpectrum(i)
        s= removespikes(s,thresholds=[3,3])
        
        s.plot()
        fluorescencedots = np.append(fluorescencedots, s.calc_area((480,620)))
    
    anth =  RamanSpectrum('160217/Raman160217/160217_25.txt')
    antharea = anth.calc_area((500,720))
    
    nliq = 1.333
    nE = 1.359
    
    fluoryields = fluorescencedots*0.65*(1-10**(-absorbance_anth))*nliq**2/nE**2 /antharea/(1-10**(-absorbances_dots)) 
    
    print 'oleate capped dots QY:', fluoryields[0]
    fluoryields_PPA = transpose(fluoryields[1:].reshape((-1,4)))
    figure()
    plot([5,7,8,9,11],fluoryields_PPA[0])
    plot([5,7,8,9,11],fluoryields_PPA[1])
    plot([5,7,8,9,11],fluoryields_PPA[2])
    plot([5,7,8,9,11],fluoryields_PPA[3])
    print transpose(filenames.reshape((-1,4)))[0]
    
    
    return fluoryields_PPA.flatten()

def Feb17CdSekr():
    timeconstants_s1 = Feb17CdSeTCSPC(2)  
    tau_ave = array(sum(timeconstants_s1['tau']**2*timeconstants_s1['amp'],axis=1)/sum(timeconstants_s1['tau']*timeconstants_s1['amp'],axis=1))
    
    
    print 'oleate ks', 1/timeconstants_s1['tau'][0]
    print 'oleate kr, knr1,knr2', 0.026/tau_ave[0], 1/timeconstants_s1['tau'][0]-0.026/tau_ave[0]
    tau_ave = tau_ave[array([1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23])]
    QYs =  Feb17CdSeQYs()
    
    kr = (QYs/tau_ave).reshape((4,5))
    print timeconstants_s1.shape, QYs.shape
    knrs = 1/timeconstants_s1['tau'][array([1,2,3,4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23])]-np.transpose([(QYs/tau_ave)]*2)
    print knrs
    kr_ave = np.mean(kr,axis = 0)
    kr_std = np.std(kr,axis = 0)
    close('all')
    figure()
    plt.errorbar([5,7,8,9,11],kr_ave, yerr=kr_std)
  
    xlabel('pH')
    ylabel('k$_r$ (ns$^{-1}$)')
    
    print 'K_nr', knrs
    
    return (kr_ave, kr_std)