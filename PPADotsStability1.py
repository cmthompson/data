# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:53:24 2015

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

universalfilelist = [ 'dot1d_N2pH7', 'dot1d_N2pH9', 'dot1d_N2pH11', 
                    'dot1d_airpH7', 'dot1d_airpH9', 'dot1d_airpH11',
                    'dot2d_N2pH7', 'dot2d_N2pH9','dot2d_N2pH11',
                    'dot2d_airpH7', 'dot2d_airpH9','dot2d_airpH11' ]
                    
O2_quench_correction = 1+0.00145*158
#def PPAcappeddotsstability():
#    
#    indexlist = [ 'dot1N2pH7', 'dot1N2pH9', 'dot1N2pH11', 
#                    'dot1aipH7', 'dot1aipH9', 'dot1aipH11',
#                    'dot2N2pH7', 'dot2N2pH9','dot2N2pH11',
#                    'dot2aipH7', 'dot2aipH9','dot2aipH11', ]
#    data0 = array([Sept6_Fluorescence_Dotsinwater()[0]])
#
#    data0 =numpy.append(data0, array([Sept7_Fluorescence_Dotsinwater()[0]]),axis=0)
#    data0 =numpy.append(data0, array([Sept8_Fluorescence_Dotsinwater()[0]]),axis=0)
#    data0 =numpy.append(data0, array([Sept9_Fluorescence_Dotsinwater()[0]]),axis=0)
#    data0 =numpy.append(data0, array([Sept10_Fluorescence_Dotsinwater()[0]]),axis=0)
#    data0 =numpy.append(data0, array([Sept11_Fluorescence_Dotsinwater()[0]]),axis=0)
#    data0 =numpy.append(data0, array([Sept13_Fluorescence_Dotsinwater()[0]]),axis=0)
#    
#    fluorescencespectra =array([Sept6_Fluorescence_Dotsinwater()[1]]) 
#   
#    fluorescencespectra =numpy.append(fluorescencespectra, array([Sept7_Fluorescence_Dotsinwater()[1]]),axis=0)
#    fluorescencespectra =numpy.append(fluorescencespectra, array([Sept8_Fluorescence_Dotsinwater()[1]]),axis=0)
#    fluorescencespectra =numpy.append(fluorescencespectra, array([Sept9_Fluorescence_Dotsinwater()[1]]),axis=0)
#    fluorescencespectra =numpy.append(fluorescencespectra, array([Sept10_Fluorescence_Dotsinwater()[1]]),axis=0)
#    fluorescencespectra =numpy.append(fluorescencespectra, array([Sept11_Fluorescence_Dotsinwater()[1]]),axis=0)
#    fluorescencespectra =numpy.append(fluorescencespectra, array([Sept13_Fluorescence_Dotsinwater()[1]]),axis=0)
#    
#    
#    fig1=figure()
#    fig1.add_subplot((131))
#    fig1.add_subplot((132))
#    fig1.add_subplot((133))
#    fluordata = pandas.DataFrame(data0,columns =indexlist ,index =[0,1,2,3,4,5,7] )
#    fluordata['dot1N2pH9'][:]=0
#
#    for i in range(3):
#        fluordata[indexlist[i]].plot(color = 'r',marker='s',ax=fig1.axes[i])
#        fluordata[indexlist[3+i]].plot(color = 'k',marker ='s',ax=fig1.axes[i])
#        fluordata[indexlist[6+i]].plot(color = 'r',marker='s',ax=fig1.axes[i])
#        fluordata[indexlist[9+i]].plot(color = 'k',marker ='s',ax=fig1.axes[i])
#       
#        fig1.axes[i].legend()
#        fig1.axes[i].set_ylim(0,0.0012)
#    numberofdays = 7
#  
#    figfluorspectra = figure(figsize=(10,10))
#
#    for r in range(6):  ### the different samples
#        subplot(320+r+1)
#        
#        title(indexlist[r])
#        for i in range(6):
#
#            try:
#                x=arange(400,680)
#                plot(x,fluorescencespectra[i,r],'-')
#            except:
#                print 'could not plot dots1', day, r
#        #ylim(0,0.2)
#        
#        legend(['day0','day1','day2','day3','day4','day5','day7'], fontsize=10)
#    figure(figsize=(10,10))
#    for r in range(6,12):  ### the different samples
#        subplot(320+r-6+1)
#        title(indexlist[r])
#        for day in range(numberofdays):
#            
#            try:
#                x=arange(400,680)
#                plot(x,fluorescencespectra[day,r],'-')
#            except:
#                print 'could not plot dots2, day', day,', spectrum', r-6
#      #  ylim(0,0.2)
#        legend(['day0','day1','day2','day3','day4','day5','day7'], fontsize=10)
#    
#    return 0
#    
#def UVVisfigure():    
#    indexlist = [ 'dot1N2pH7', 'dot1N2pH9', 'dot1N2pH11', 
#                    'dot1aipH7', 'dot1aipH9', 'dot1aipH11',
#                    'dot2N2pH7', 'dot2N2pH9','dot2N2pH11',
#                    'dot2aipH7', 'dot2aipH9','dot2aipH11', ]
#    
#    
#    UVVisspectra = array([Sept6_UVVis_Dotsinwater()[1]])
#        
#  ## data0 array goes [day, spectrum, frequency]
#    
#    UVVisspectra =numpy.append(UVVisspectra, array([Sept7_UVVis_Dotsinwater()[1]]),axis=0)
#    UVVisspectra =numpy.append(UVVisspectra, array([Sept8_UVVis_Dotsinwater()[1][:,100:]]),axis=0)
#    UVVisspectra =numpy.append(UVVisspectra, array([Sept9_UVVis_Dotsinwater()[1]]),axis=0)
#    UVVisspectra =numpy.append(UVVisspectra, array([Sept10_UVVis_Dotsinwater()[1]]),axis=0)
#    UVVisspectra =numpy.append(UVVisspectra, array([Sept11_UVVis_Dotsinwater()[1]]),axis=0)
#    UVVisspectra =numpy.append(UVVisspectra, array([Sept13_UVVis_Dotsinwater()[1]]),axis=0)
#    UVVisspectra-=reshape(UVVisspectra[:,:,0],(7,13,1))
#    
#    
#    numberofdays = UVVisspectra.shape[0]
#    fig2 = figure(figsize=(10,10))
#    for r in range(6):  ### the different samples
#        subplot(320+r+1)
#        title(indexlist[r])
#        for i in range(numberofdays):
#            plot(arange(500,199,-1),UVVisspectra[i,r],'-',label='day'+str(i))
#        ylim(0,0.1)
#        legend(loc=2, fontsize = 10)
#        xlim(360,480)
#    savefig('/home/chris/Dropbox/DataWeiss/150913/UVspectradots1.jpg')
#    figure(figsize=(10,10))
#    for r in range(6,12):  ### the different samples
#        subplot(320+r-6+1)
#        title(indexlist[r])
#        for day in range(numberofdays):
#            plot(arange(500,199,-1),UVVisspectra[day,r],'-',label='day'+str(day))
#        ylim(0,0.1)
#        xlim(360,480)
#        legend(loc=2)
#    savefig('/home/chris/Dropbox/DataWeiss/150913/UVspectradots2.jpg')
#    
#    fig2 = figure(figsize=(10,10))
#    for r in range(6):  ### the different samples
#        subplot(320+r+1)
#        title(indexlist[r])
#        for i in range(numberofdays):
#            if i==6:
#                daylabel=7
#            else:
#                daylabel=i
#            plot(arange(500,199,-1),UVVisspectra[i,r],'-',label='day'+str(daylabel))
#        ylim(0,1)
#        legend(loc=2, fontsize = 10)
#        
#    savefig('/home/chris/Dropbox/DataWeiss/150913/UVspectradots1showDMF.jpg')
#    figure(figsize=(10,10))
#    for r in range(6,12):  ### the different samples
#        subplot(320+r-6+1)
#        title(indexlist[r])
#        for day in range(numberofdays):
#            if day==6:
#                daylabel=7
#            else:
#                daylabel=day
#            plot(arange(500,199,-1),UVVisspectra[day,r],'-',label='day'+str(daylabel))
#        ylim(0,1)
#        legend(loc=2,fontsize = 10)
#    savefig('/home/chris/Dropbox/DataWeiss/150913/UVspectradots2shwoDMF.jpg')
#    
#    return 0
#
#
#
#    
#def Sept6_UVVis_Dotsinwater():
#
#    a = loadtxt('/home/chris/Dropbox/DataWeiss/150906/PPA dots in water for stability study.csv', delimiter = ',',usecols=[0,4,5,6,9,10,11,14,15,16,19,20,21,23], unpack = True, skiprows = 1)
#    UVVis = ndarray((0,len(a[0])))    
#    f = open('/home/chris/Dropbox/DataWeiss/150906/PPA dots in water for stability study.csv','rb')   
#    z = f.readline()
#    f.close()
#    z= z.split(',')
#    if False:  ########### This lists the samples
#        for i in [0,4,5,6,9,10,11,14,15,16,19,20,21]:
#            print z[i]
#        for i in a[1:]:
#            i[:]-=i[0]
#            plot(a[0],i)
#
#    ######## THE FILES LOOK MIXED UP BECAUSE I CHNAGED THE "O2" SAMPLES TO BE "N2" 
#    UVVis=a[1:,:]-transpose([numpy.min(a[1:],axis = 1)]) ### bring minimum to zero
#    x = argmin(abs(a[0]-350))
#    absorbancevalues = a[1:13,x]  
#    anthraceneabsorbance = a[-1,x]
#    
#  
#    return (absorbancevalues,UVVis,anthraceneabsorbance)
def correctionfactorforunfilledcuvettes():
#     print 'The data from september6 was incorrect, since the cuvettes were not completely full, so that the fluorescence from the'
#     print 'samples was too high.  Below are the points collected for Sept7 in the "high" position (so that liquid covered the whole probe area"'
#     print 'and "low" such that the sample was in the same position as Sept6.'
#     print 'the ratio of "hi"to"low" is'
     os.chdir('/home/chris/Dropbox/DataWeiss/150907/150907fluor')
    
     
     filelist = ['dot1d1N2pH11h9', 'dot1d1N2pH11l8', 
                    'dot1d1N2pH7hi1','dot1d1N2pH7lo1',
                    'dot1d1N2pH7hi4', 'dot1d1N2pH7lo2',
                    'dot1d1N2pH9hi7', 'dot1d1N2pH9lo6',
                    'dot1d1aipH11h6', 'dot1d1aipH11l5',
                    'dot1d1aipH7h1', 'dot1d1aipH7l1', 
                    'dot1d1aipH9h4', 'dot1d1aipH9l3', 
                    'dot2d1N2pH11h5', 'dot2d1N2pH11l4',
                    'dot2d1N2pH7h1', 'dot2d1N2pH7l10', 
                    'dot2d1N2pH9h1', 'dot2d1N2pH9l1', 
                    'dot2d1aipH11h6', 'dot2d1aipH11l5',
                    'dot2d1aipH7h1', 'dot2d1aipH7l1',
                    'dot2d1aipH9h4', 'dot2d1aipH9l1']
     
     thelist = list()
     for i in range(len(filelist)/2):
         if filelist[2*i][0:11]!=filelist[2*i+1][0:11]:
             print filelist[2*i],filelist[2*i+1]
             continue
         else:
             a = loadtxt(filelist[2*i],delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
             hi = RamanSpectrum(pandas.Series(a[1],a[0]))
             b = loadtxt(filelist[2*i+1],delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
             lo = RamanSpectrum(pandas.Series(b[1],b[0]))
             hi.plot()
             lo.plot()
             hidivlo= hi.calc_area((410,473))/lo.calc_area((410,473))
             thelist.append(hidivlo)
     thelist=array(thelist)
#     print mean(thelist)
#     print 'with a standard deviation of ', std(thelist)
#     print  "The values are returned here"
     ## DATA FROM SEPT6 SHOULD BE MULTIPLIED BY THIS NUMBER TO GET THE CORRECT VALUE
     return numpy.mean(thelist)
#     
#     
#def Sept6_Fluorescence_Dotsinwater():
#     folder = '/home/chris/Dropbox/DataWeiss/150906/150906fluor/'
#    
#     os.chdir(folder)
#     a = loadtxt('aminoanthracene.dat',delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
#     x = argmin(abs(a[0]-450))  ###Normalizing to value at 450 nm
#     Ianth = a[1,x]-a[1,0]
#     
#     print 'anthracene standard for'+folder.split('/')[-1]+' at 450 nm:', Ianth
#
#     (absvalues, UVVis,Aanth) =Sept6_UVVis_Dotsinwater()
#     
#     correctionfactor = 0.9034907#correctionfactorforunfilledcuvettes()
#     print 'using correction factor for the unfilled cuvettes on Sept6 data:', correctionfactor
#    
#         
#    ######## THE FILES LOOK MIXED UP BECAUSE I CHNAGED THE "O2" SAMPLES TO BE "N2"                
#     filelist = ['dot1_O2_pH7.dat','dot1_O2_pH9.dat','dot1_O2_pH11.dat',
#                'dot1_N2_pH7.dat','dot1_N2_pH9.dat','dot1_N2_pH11.dat',
#                'dot2_O2_pH7.dat','dot2_O2_pH9.dat','dot2_O2_pH11.dat',
#                'dot2_N2_pH7.dat','dot2_N2_pH9.dat','dot2_N2_pH11.dat',]
#  
#     
#     arealist = list()
#     fluorescence = ndarray((0,280))
#     for i in range(len(filelist)):
#         
#         if os.path.exists(filelist[i]):
#             a = loadtxt(filelist[i],delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
#             a[1]-=min(a[1])
#             a[1]*=0.000380*(1-10**(-Aanth))/(1-10**(-absvalues[i]))/Ianth
#             hi = RamanSpectrum(pandas.Series(a[1],a[0]))
#             hi.truncate(after = 681)
#              
#             arealist.append(hi.calc_area((410,473)))
#             fluorescence = append(fluorescence ,array([a[1]]), axis = 0)
#         else:
#             print 'problem with',folder,filelist[i]
#             arealist.append(0)
#             fluorescence = append(fluorescence,zeros((1,280)),axis=0)
#     
#     areas = array(arealist)
#     for i in fluorescence:
#         plot(i)
#     legend(list(str(i) for i in [0,1,2,3,4,5,7]))
#     return (areas,fluorescence)
#     
#
#
#def Sept7_UVVis_Dotsinwater():
#
#    return workupUVVis('/home/chris/Dropbox/DataWeiss/150907/150907UVVis.csv','p')
#
#     
#def Sept7_Fluorescence_Dotsinwater():
#    #########2-aminoanthracene normalization
#     folder = '/home/chris/Dropbox/DataWeiss/150907/150907fluor/'
#     os.chdir(folder)
#     a = loadtxt('aminoanthracene.dat',delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
#     x = argmin(abs(a[0]-450))  ###Normalizing to value at 450 nm
#     Ianth = a[1,x]-a[1,0]
#     
#     print 'anthracene standard for'+folder.split('/')[-1]+' at 450 nm:', Ianth
#
#     (absvalues, UVVis,Aanth) =Sept7_UVVis_Dotsinwater()
#     
#     
#     filelist = [ 'dot1d1N2pH7hi1', 'dot1d1N2pH9hi7', 'dot1d1N2pH11h9', 
#                    'dot1d1aipH7h1', 'dot1d1aipH9h4', 'dot1d1aipH11h6',
#                    'dot2d1N2pH7h1', 'dot2d1N2pH9h1','dot2d1N2pH11h5',
#                    'dot2d1aipH7h1', 'dot2d1aipH9h4','dot2d1aipH11h6', ]
#     arealist = list()
#     fluorescence = ndarray((12,280))
#     for i in range(len(filelist)):
#         
#         if os.path.exists(filelist[i]):
#             a = loadtxt(filelist[i],delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
#             a[1]-=min(a[1])
#             a[1]*=0.000380*(1-10**(-Aanth))/(1-10**(-absvalues[i]))/Ianth
#             hi = RamanSpectrum(pandas.Series(a[1],a[0]))
#             hi.truncate(before=400,after = 680)
#              
#             arealist.append(hi.calc_area((410,473)))
#             fluorescence[i,0:min(280,len(a[1]))] = a[1,:min(280,len(a[1]))]
#             #fluorescence = append(fluorescence ,array([a[1]]), axis = 0)
#         else:
#             print 'problem with',folder,filelist[i]
#             arealist.append(0)
#             
#     areas = array(arealist)
#     for i in fluorescence:
#         plot(i)
#     legend(list(str(i) for i in [0,1,2,3,4,5,7]))
#     return (areas,fluorescence)
#    
#         
#   
#     
#def Sept8_UVVis_Dotsinwater():
#    return workupUVVis('/home/chris/Dropbox/DataWeiss/150908/PPAstability.csv','o')
##    a = loadtxt(filename, delimiter = ',', unpack = True, skiprows = 101, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,14))
##    UVVis=a[1:,:]-transpose([numpy.min(a[1:],axis = 1)]) ### bring minimum to zero
##    x = argmin(abs(a[0]-350))
##    absorbancevalues = a[1:13,x]  
##    anthraceneabsorbance = a[-1,x]
##
##    return (absorbancevalues,UVVis,anthraceneabsorbance)
#    
#def Sept8_Fluorescence_Dotsinwater():
#     return workupfluorescence('/home/chris/Dropbox/DataWeiss/150908/150908fluor/',Sept8_UVVis_Dotsinwater,'2')
#     
#def Sept9_UVVis_Dotsinwater(_plot=False):
#     return workupUVVis('/home/chris/Dropbox/DataWeiss/150909/150909UVVis.csv','o')
##    a = loadtxt('/home/chris/Dropbox/DataWeiss/150909/150909UVVis.csv', delimiter = ',', unpack = True, skiprows = 1)
##    a = a[:13]
##    UVVis = ndarray((0,len(a[0])))    
##    f = open('/home/chris/Dropbox/DataWeiss/150909/150909UVVis.csv','rb')   
##    
##    z = f.readline()
##    f.close()
##    z= z.split(',')
##    if _plot == True:  ########### This lists the samples
##        for i in range(1,13):
##            print z[i]
##        
##            a[i][:]-=a[i][0]
##            plot(a[0],a[i],label = z[i])
##        legend()
##        ylim(0,0.5)
##    for i in a[1:]:
##        i-=i[0]
##        UVVis = append(UVVis, array([i]),axis = 0)
##    x = argmin(abs(a[0]-350))
##
##    absorbancevalues = a[1:,x]
##    return (absorbancevalues,UVVis)
#    
#
#     
#     
#def workupUVVis(filename,anthracenecolumn,_plot=False):
#    alphabet = 'abcdefghijklmnopqrstuvwxyz'
#    
#    a = loadtxt(filename, delimiter = ',', unpack = True, skiprows = 1, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,alphabet.find(anthracenecolumn)))
#    
#    
#    UVVis=a[1:,:]-transpose([a[1:,0]]) ### bring minimum to zero
#    x = argmin(abs(a[0]-350))
#    absorbancevalues = a[1:13,x]  
#    anthraceneabsorbance = a[-1,x]
#
#    if _plot:  ########### This lists the samples
#        f = open(filename,'rb')   
#        z = f.readline()
#        f.close()
#        z= z.split(',')
#        for i in range(1,13):
#            
#            a[i][:]-=a[i][0]
#            plot(a[0],a[i],label = z[i])
#        legend()
#        ylim(0,0.5)
#    
#    return (absorbancevalues,UVVis,anthraceneabsorbance)
#    
#def aminoanthracenespectra():
#     a = loadtxt('/home/chris/Dropbox/DataWeiss/150906/PPA dots in water for stability study.csv', delimiter = ',',usecols=(0,23), unpack = True, skiprows = 1)
#     b = loadtxt('/home/chris/Dropbox/DataWeiss/150907/150907UVVis.csv',delimiter=',', unpack = True,skiprows=1,usecols=[0,15])
#     c = loadtxt('/home/chris/Dropbox/DataWeiss/150908/PPAstability.csv', delimiter = ',', unpack = True, skiprows = 1, usecols=(0,14))
#     
#    
#     d = loadtxt('/home/chris/Dropbox/DataWeiss/150909/150909UVVis.csv', delimiter = ',', unpack = True, skiprows = 1,usecols=(0,14))
#     e = loadtxt('/home/chris/Dropbox/DataWeiss/150910/150910UVVis.csv', delimiter = ',', unpack = True, skiprows = 1,usecols=(0,13))
#     
#     f = loadtxt('/home/chris/Dropbox/DataWeiss/150911/150911UVVis.csv',delimiter=',', unpack = True,skiprows=1,usecols=[0,14])
#     g = loadtxt('/home/chris/Dropbox/DataWeiss/150913/150913UVVis.csv',delimiter=',', unpack = True,skiprows=1,usecols=[0,14])
#         
#     for i in (a,b,c,d,e,f,g):
#         i[1]-=i[1,0]
#         plot(i[0],i[1])
#     legend(['0','1','2','3','4','5','7'])
#     
#     return 0
#     
#def workupfluorescence(folder,UVVisfunction,day):
#    #########2-aminoanthracene normalization
#     os.chdir(folder)
#     a = loadtxt('aminoanthracene.dat',delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
#     x = argmin(abs(a[0]-450))  ###Normalizing to value at 450 nm
#     Ianth = a[1,x]-a[1,0]
#     
#     print 'anthracene standard for'+folder.split('/')[-1]+' at 450 nm:', Ianth
#
#     (absvalues, UVVis,Aanth) = UVVisfunction()
#     filelist = list((x.replace('_', day) for x in universalfilelist))
#  
#     
#     arealist = list()
#     fluorescence = ndarray((0,280))
#     for i in range(len(filelist)):
#         
#         if os.path.exists(filelist[i]+'.dat'):
#             a = loadtxt(filelist[i]+'.dat',delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
#             a[1]-=min(a[1])
#             a[1]*=0.000380*(1-10**(-Aanth))/(1-10**(-absvalues[i]))/Ianth
#             hi = RamanSpectrum(pandas.Series(a[1],a[0]))
#             hi.truncate(after = 681)
#              
#             arealist.append(hi.calc_area((410,473)))
#             fluorescence = append(fluorescence ,array([a[1]]), axis = 0)
#         else:
#             print 'problem with',folder,filelist[i]
#             arealist.append(nan)
#             fluorescence = append(fluorescence,zeros((1,280)),axis=0)
#     
#     areas = array(arealist)
#    
#     return (areas,fluorescence)
#     
#def Sept9_Fluorescence_Dotsinwater():
#    #########2-aminoanthracene normalization
#     return workupfluorescence('/home/chris/Dropbox/DataWeiss/150909/150909fluor',Sept9_UVVis_Dotsinwater,'3')
#     
#     
#def Sept10_UVVis_Dotsinwater(_plot=False):
#    return workupUVVis('/home/chris/Dropbox/DataWeiss/150910/150910UVVis.csv','n')
#    
#    
#def Sept10_Fluorescence_Dotsinwater():
#    #########2-aminoanthracene normalization
#     return workupfluorescence('/home/chris/Dropbox/DataWeiss/150910/150910fluor', Sept10_UVVis_Dotsinwater,'4')
#     
#    
#def Sept11_UVVis_Dotsinwater(_plot=False):
#    return workupUVVis('/home/chris/Dropbox/DataWeiss/150911/150911UVVis.csv','m',_plot = _plot)
#    
#    
#def Sept11_Fluorescence_Dotsinwater():
#    #########2-aminoanthracene normalization
#     return workupfluorescence('/home/chris/Dropbox/DataWeiss/150911/150911fluor',Sept11_UVVis_Dotsinwater,'5')
#
#def Sept13_UVVis_Dotsinwater(_plot=False):
#    return workupUVVis('/home/chris/Dropbox/DataWeiss/150913/150913UVVis.csv','m',_plot = _plot)
#    
#    
#def Sept13_Fluorescence_Dotsinwater():
#    #########2-aminoanthracene normalization
#     return workupfluorescence('/home/chris/Dropbox/DataWeiss/150913/150913fluor',Sept13_UVVis_Dotsinwater,'7')

def aminoanthracene_fluorescenceyield_determination():
    amino = loadtxt('/home/chris/Dropbox/DataWeiss/150929/150929fluor/aminoanthracene.dat',delimiter = '\t', unpack = True, usecols=(0,3),skiprows = 1)
    anthracene = loadtxt('/home/chris/Dropbox/DataWeiss/150929/150929fluor/anthracene.dat',delimiter = '\t', unpack = True, usecols=(0,3),skiprows = 1)
    amino[1]-=min(amino[1])
    anthracene[1]-=min(anthracene[1])   
    x = argmin(abs(450-amino[0]))
    Iamino = amino[1,x]
    
    x = argmin(abs(420-anthracene[0]))
    Ianth_420 = anthracene[1,x]
    
    plot(amino[0],amino[1])
    plot(anthracene[0],anthracene[1])

    uv = loadtxt('/home/chris/Dropbox/DataWeiss/150929/150929UVVis.csv', skiprows = 1, unpack=True, delimiter = ',',usecols=(0,7,8))
    Absspecanth = RamanSpectrum(pandas.Series(uv[2][::-1],uv[0][::-1]))
    Absspecamino = RamanSpectrum(pandas.Series(uv[1][::-1],uv[0][::-1]))
    
    Absspecamino-=min(Absspecamino)
    Absspecanth-=min(Absspecanth)
    Aanth = Absspecanth[350]
    Aamino = Absspecamino[350]
    
    O2quench_correction = (1+0.00145*158)
    
    QY = 0.27/O2quench_correction*Iamino*(1-10**(-Aanth))/(Ianth_420*76.9811)/(1-10**(-Aamino))
    print 'dQ/dlambda for aminoanthracene at 450 is',  QY
    figure()
    ax=subplot(111)
    r = indivQY('/home/chris/Dropbox/DataWeiss/150929/150929UVVis.csv','h','i', '/home/chris/Dropbox/DataWeiss/150929/150929fluor/aminoanthracene.dat','/home/chris/Dropbox/DataWeiss/150929/150929fluor/anthracene.dat',fluorescencerange = (450,455), UVVisplot=ax, fluorplot = ax, output='fluorescencespectrum' )
    print r[450]    
    return 0
    
def indivQY(UVVisfile, UVViscolumn, anthracenecolumn, fluorescencefile, anthracenefluorescencefile,
            UVVisplot=None, fluorplot=None,
            fluorescencerange = (410,473),
            nliq=1.33,
            day=0, label=None,color = 'k', output = None):
    print '-------------------------------------'
    print 'calculating fluorescence yield for', label,'file', fluorescencefile
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    if len(UVViscolumn)==1:
        numuvviscolumn = alphabet.find(UVViscolumn)
    elif len(UVViscolumn)==2:
        numuvviscolumn = alphabet.find(UVViscolumn[0])*26+alphabet.find(UVViscolumn[1])
    
        
    
    a = loadtxt(UVVisfile, delimiter = ',', unpack = True, skiprows = 1,usecols=(0,numuvviscolumn,alphabet.find(anthracenecolumn)))
    a[1:]-=transpose([a[1:,0]])
    

    anthracene=RamanSpectrum(pandas.Series(a[2][::-1],a[0][::-1]))
    dot=RamanSpectrum(pandas.Series(a[1][::-1],a[0][::-1]))
    
    anthraceneabsorbance350= anthracene[350]#(anthracene[374]-anthracene[389])*0.6735# 
    absvalues = dot[350]
    
    nE = 1.5
    nQ = 1.44
    nW = 1.33  ## refractive index water    
    
    TinEtOH= 1 - ((nQ-nE)/(nQ+nE))**2 
    ToutEtOH= 1 - ((nE-nQ)/(nE+nQ))**2 
    TinH2O =  1 - ((nQ-nliq)/(nQ+nliq))**2 
    ToutH2O = 1 - ((nW-nliq)/(nW+nliq))**2 
    
    
    a = loadtxt(anthracenefluorescencefile,delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
    a[1]-=a[1,-1]
    anthracenefluorescence=RamanSpectrum(pandas.Series(a[1],a[0]))
    
     ###Normalizing to value of anthracene at 420 nm The area for the anthracene fluorescence is related to this value by 78.203    
    anthracenefluorescencearea = anthracenefluorescence[450]#*78.2032212661
    print 'anthracene fluorescence area=', '%.2E' % anthracenefluorescencearea
   
    oneminusTdot = 1-10**(-absvalues)   ##### gives the fraction of photons absorbed by dots
    
    oneminusT_anthracene350 =1-10**(-anthraceneabsorbance350)
    print 'anthracene absorbance at 350 nm:', anthraceneabsorbance350,'. Fraction photons absorbed:', oneminusT_anthracene350
    print 'dot absorbance at 350 nm:', absvalues, '. Fraction photons absorbed:', oneminusTdot

    
    a = loadtxt(fluorescencefile,delimiter='\t', unpack = True,skiprows=1, usecols = (0,3))
    hi = RamanSpectrum(pandas.Series(a[1],a[0]))
    hi[:]-=hi[407]
    hi[:]*=0.0001484*oneminusT_anthracene350*nliq**2/nE**2 /anthracenefluorescencearea/ oneminusTdot 
    
    
    dotfluorescencearea = hi.calc_area(fluorescencerange,fill=False)
    
    
    ## quantum yield of dots using 0.27 as QY for anthracene with o2 quenching corrrection
    print  'fluorescence (bande edg) yield of dot', dotfluorescencearea
    
    if UVVisplot is not None:
       # anthracene.plot(ax=UVVisplot)#plot(a[0],anthracene)
        dot.plot(ax=UVVisplot,label=label)
        UVVisplot.plot(350,dot[350],'s')
    if fluorplot is not None:
        hi.plot(ax = fluorplot,label=label)
        
        #anthracenefluorescence.plot(ax = fluorplot,label=label)
    if output=='fluorescencespectrum':
        return hi
    elif output=='uvvisspectrum':
        return dot
    else:
        return  dotfluorescencearea
    

def trial1_Sept6through13():
    ## day 0
    sample1=list()
    sample2=list()
    sample3=list()
    sample4=list()
    sample5=list()
    sample6=list()
    sample7=list()
    sample8=list()
    sample9=list()
    sample10=list()
    sample11=list()
    sample12=list()
    
    figure()
    Aax1 = subplot(231)
    Aax2 = subplot(232)
    Aax3 = subplot(233)
    Aax4 = subplot(234)
    Aax5 = subplot(235)
    Aax6 = subplot(236)
    
    figure()
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(234)
    ax5 = subplot(235)
    ax6 = subplot(236)
    ##day0
  
    uvvisfile = '/home/chris/Dropbox/DataWeiss/150906/PPA dots in water for stability study.csv'  #, delimiter = ',',usecols=[0,4,5,6,9,10,11,14,15,16,19,20,21,23]
    anthracenefile = '/home/chris/Dropbox/DataWeiss/150906/150906fluor/aminoanthracene.dat'
    fluorfolder = '/home/chris/Dropbox/DataWeiss/150906/150906fluor'
    sample1.append(0.907*indivQY(uvvisfile, 'e','i',fluorfolder+'/dot1_O2_pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0.907*indivQY(uvvisfile, 'f','i',fluorfolder +'/dot1_O2_pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(0.907*indivQY(uvvisfile, 'g','i',fluorfolder +'/dot1_O2_pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(0.907*indivQY(uvvisfile, 'j','i',fluorfolder +'/dot1_N2_pH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(0.907*indivQY(uvvisfile, 'k','i',fluorfolder +'/dot1_N2_pH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(0.907*indivQY(uvvisfile, 'l','i',fluorfolder +'/dot1_N2_pH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(0.907*indivQY(uvvisfile, 'o','i',fluorfolder+'/dot2_O2_pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(0.907*indivQY(uvvisfile, 'p','i',fluorfolder +'/dot2_O2_pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(0.907*indivQY(uvvisfile, 'q','i',fluorfolder +'/dot2_O2_pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(0.907*indivQY(uvvisfile, 't','i',fluorfolder +'/dot2_N2_pH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(0.907*indivQY(uvvisfile, 'u','i',fluorfolder +'/dot2_N2_pH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(0.907*indivQY(uvvisfile, 'v','i',fluorfolder +'/dot2_N2_pH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    
    
    ##day1
    uvvisfile = '/home/chris/Dropbox/DataWeiss/150907/150907UVVis.csv'  #, delimiter = ',',usecols=[0,4,5,6,9,10,11,14,15,16,19,20,21,23]
    anthracenefile = '/home/chris/Dropbox/DataWeiss/150907/150907fluor/aminoanthracene.dat'
    fluorfolder = '/home/chris/Dropbox/DataWeiss/150907/150907fluor/'
    sample1.append(indivQY(uvvisfile, 'b','p',fluorfolder+'dot1d1N2pH7hi1',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
   # sample2.append(indivQY(uvvisfile, 'c','p',fluorfolder + 'dot1d1N2pH9hi7',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','p',fluorfolder +'dot1d1N2pH11h9', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','p',fluorfolder +'dot1d1aipH7h1', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','p',fluorfolder +'dot1d1aipH9h4', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','p',fluorfolder +'dot1d1aipH11h6',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','p',fluorfolder+'dot2d1N2pH7h1',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','p',fluorfolder + 'dot2d1N2pH9h1',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','p',fluorfolder +'dot2d1N2pH11h5', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','p',fluorfolder +'dot2d1aipH7h1', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','p',fluorfolder +'dot2d1aipH9h4', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','p',fluorfolder +'dot2d1aipH11h6',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
  
    ##day2

     
    uvvisfile = '/home/chris/Dropbox/DataWeiss/150908/PPAstability.csv'
    anthracenefile = '/home/chris/Dropbox/DataWeiss/150908/150908fluor/aminoanthracene.dat'
    fluorfolder = '/home/chris/Dropbox/DataWeiss/150908/150908fluor/'
    sample1.append(indivQY(uvvisfile, 'b','o',fluorfolder+'dot1d2N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0)#ivQY(uvvisfile, 'c','o',fluorfolder + 'dot1d2N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','o',fluorfolder +'dot1d2N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','o',fluorfolder +'dot1d2airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','o',fluorfolder +'dot1d2airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','o',fluorfolder +'dot1d2airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','o',fluorfolder+'dot2d2N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','o',fluorfolder + 'dot2d2N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','o',fluorfolder +'dot2d2N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','o',fluorfolder +'dot2d2airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','o',fluorfolder +'dot2d2airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','o',fluorfolder +'dot2d2airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
   
    ##day3
    uvvisfile = '/home/chris/Dropbox/DataWeiss/150909/150909UVVis.csv'
    anthracenefile = '/home/chris/Dropbox/DataWeiss/150909/150909fluor/aminoanthracene.dat'
    fluorfolder = '/home/chris/Dropbox/DataWeiss/150909/150909fluor/'
    sample1.append(indivQY(uvvisfile, 'b','o',fluorfolder+'dot1d3N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
   # sample2.append(0)#ivQY(uvvisfile, 'c','o',fluorfolder + 'dot1d2N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','o',fluorfolder +'dot1d3N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','o',fluorfolder +'dot1d3airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','o',fluorfolder +'dot1d3airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','o',fluorfolder +'dot1d3airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','o',fluorfolder+'dot2d3N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','o',fluorfolder + 'dot2d3N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','o',fluorfolder + 'dot2d3N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','o',fluorfolder +'dot2d3airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','o',fluorfolder +'dot2d3airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','o',fluorfolder +'dot2d3airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
#    ##day4
 
    uvvisfile ='/home/chris/Dropbox/DataWeiss/150910/150910UVVis.csv'
    anthracenefile = '/home/chris/Dropbox/DataWeiss/150910/150910fluor/aminoanthracene.dat'
    fluorfolder = '/home/chris/Dropbox/DataWeiss/150910/150910fluor/'
    sample1.append(indivQY(uvvisfile, 'b','n',fluorfolder+'dot1d4N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
   # sample2.append(0)#ivQY(uvvisfile, 'c','n',fluorfolder + 'dot1d2N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','n',fluorfolder +'dot1d4N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','n',fluorfolder +'dot1d4airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','n',fluorfolder +'dot1d4airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','n',fluorfolder +'dot1d4airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','n',fluorfolder+'dot2d4N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','n',fluorfolder +'dot2d4N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','n',fluorfolder +'dot2d4N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','n',fluorfolder +'dot2d4airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','n',fluorfolder +'dot2d4airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','n',fluorfolder +'dot2d4airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
#    ##day5
#    day=5
    uvvisfile ='/home/chris/Dropbox/DataWeiss/150911/150911UVVis.csv'
    anthracenefile = '/home/chris/Dropbox/DataWeiss/150911/150911fluor/aminoanthracene.dat'
    fluorfolder = '/home/chris/Dropbox/DataWeiss/150911/150911fluor/'
    sample1.append(indivQY(uvvisfile, 'b','m',fluorfolder+'dot1d5N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0)#ivQY(uvvisfile, 'c','m',fluorfolder + 'dot1d2N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','m',fluorfolder +'dot1d5N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','m',fluorfolder +'dot1d5airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','m',fluorfolder +'dot1d5airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','m',fluorfolder +'dot1d5airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','m',fluorfolder+'dot2d5N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','m',fluorfolder + 'dot2d5N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','m',fluorfolder +'dot2d5N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','m',fluorfolder +'dot2d5airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','m',fluorfolder +'dot2d5airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','m',fluorfolder +'dot2d5airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
#     
#      ##day7
    day=7
    uvvisfile ='/home/chris/Dropbox/DataWeiss/150913/150913UVVis.csv'
    anthracenefile = '/home/chris/Dropbox/DataWeiss/150913/150913fluor/aminoanthracene.dat'
    fluorfolder = '/home/chris/Dropbox/DataWeiss/150913/150913fluor/'
    sample1.append(0)#indivQY(uvvisfile, 'b','m',fluorfolder+'dot1d7N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    #sample2.append(0)#ivQY(uvvisfile, 'c','m',fluorfolder + 'dot1d2N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample3.append(indivQY(uvvisfile, 'd','m',fluorfolder +'dot1d7N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample4.append(indivQY(uvvisfile, 'e','m',fluorfolder +'dot1d7airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample5.append(indivQY(uvvisfile, 'f','m',fluorfolder +'dot1d7airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample6.append(indivQY(uvvisfile, 'g','m',fluorfolder +'dot1d7airpH11.dat', anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
    sample7.append(indivQY(uvvisfile, 'h','m',fluorfolder+'dot2d7N2pH7.dat',  anthracenefile,UVVisplot = Aax1,fluorplot = ax1,day = 0,label = 'N2pH7',color = 'r'))
    sample8.append(indivQY(uvvisfile, 'i','m',fluorfolder + 'dot2d7N2pH9.dat',  anthracenefile,UVVisplot = Aax2,fluorplot = ax2,day = 0,label = 'N2pH9',color = 'r'))
    sample9.append(indivQY(uvvisfile, 'j','m',fluorfolder +'dot2d7N2pH11.dat', anthracenefile,UVVisplot = Aax3,fluorplot = ax3,day = 0,label = 'N2pH11',color = 'r'))
    sample10.append(indivQY(uvvisfile, 'k','m',fluorfolder +'dot2d7airpH7.dat', anthracenefile,UVVisplot = Aax4,fluorplot = ax4,day = 0,label = 'airpH7',color = 'r'))
    sample11.append(indivQY(uvvisfile, 'l','m',fluorfolder +'dot2d7airpH9.dat', anthracenefile,UVVisplot = Aax5,fluorplot = ax5,day = 0,label = 'airpH9',color = 'r'))
    sample12.append(indivQY(uvvisfile, 'm','m',fluorfolder +'dot2d7airpH11.dat',anthracenefile,UVVisplot = Aax6,fluorplot = ax6,day = 0,label = 'airpH11',color = 'r'))
 
    
    figure()
    ax1 = subplot(131)
    ax2 = subplot(132)
    ax3=subplot(133)
    days = [0,1,2,3,4,5,7]
    ax1.plot(days,sample1,'rs-',label = 'N2pH7')
   # plot(days,sample2,'rs-',label = 'N2pH9')
    ax3.plot(days,sample3,'rs-',label = 'N2pH11')
    ax1.plot(days,sample4,'ks-',label = 'airpH7')
    ax2.plot(days,sample5,'ks-',label = 'airpH9')
    ax3.plot(days,sample6,'ks-',label = 'airpH11')
    ax1.plot(days,sample7,'ro-',label = 'N2pH7')
    ax2.plot(days,sample8,'ro-',label = 'N2pH9')
    ax3.plot(days,sample9,'ro-',label = 'N2pH11')
    ax1.plot(days,sample10,'ko-',label = 'airpH7')
    ax2.plot(days,sample11,'ko-',label = 'airpH9')
    ax3.plot(days,sample12,'ko-',label = 'airpH11')
    ax1.set_ylim(0,0.00075)
    ax2.set_ylim(0,0.00075)
    ax3.set_ylim(0,0.00075)
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax1.set_ylabel('band edge QY')
    ax1.set_xlabel('day')
    ax2.set_xlabel('day')
    ax3.set_xlabel('day')
    return 0    