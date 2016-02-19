# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 09:02:01 2016

@author: chris
"""

from ramanTools.OrcaTools import NWChemDOS, NWchemTDDFT
import numpy as np
from ramanTools.RamanSpectrum import *
from ramanTools.OrcaTools import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import numpy as np

import pickle
from makeacluster import pointsonsphere

import load_cube
from math import factorial as fac
#from UVVistools import *
import weissdatavariables

def DFTresultschargedclusters(newonly = False):
    """workup ground  state DFT data from computations on Cd29S29 cluster with different numbers of point charges
    on it.  You do not need to use this function every time - it will save the data.  To view, use next function."""
    os.chdir('/home/chris/Dropbox/DataWeiss/160120/chargeresults')
    i = os.listdir('.')
    f=figure()
    s=0
   
    data = np.recarray((0,), dtype=[('name', '|S20'),
                                   ('charge',float),
                                    ('numcharges',int),
                                    ('distance',float),
                                    ('energy', float),
                                    ('gap', float), 
                                    ('homo', float),
                                    ('lumo', float),
                                    ('lumoaveD', float),
                                    ('homoaveD', float),
                                    ('homolumooverlap',float)])
    
   
    
    for folder in [os.path.basename(x[0]) for x in os.walk('.')][1:]:
        print folder
        if newonly ==True and folder not in os.listdir():
            continue
        else:    
            try:
                q = int(folder[folder.index('Q')+1:folder.index('d')])
            except ValueError:
                continue
            numcharges = int(folder[folder.index('d')+1:folder.index('r')])
            distance = float(folder[folder.index('r')+1:folder.index('_')])
    
            ax = None#f.add_subplot(330+s)
            if folder+'.log' not in os.listdir(folder):
                print "error.  log file not found"
            else:
                a = NWChemDOS(folder+'/'+folder+'.log',axis = ax, _plot=False)
                if 'lumo.cube' in os.listdir(folder):
                    lumocube= load_cube.CUBE(folder+'/lumo.cube')
                    lumoave = lumocube.average_distance()
                else: 
                    print "error: lumo file not found in", folder
                    lumoave = 0
                if 'homo.cube' in os.listdir(folder):
                
                    homocube= load_cube.CUBE(folder+'/homo.cube')
                    homoave = homocube.average_distance()
                    homolumoov = homocube.overlap(lumocube)# sum(homocube.data**2*lumocube.data**2)/sqrt(sum(homocube.data**4)*sum(lumocube.data**4))
                else:
                    print "error: homo file not found in", folder
                    homoave=0
            
                s+=1
               
                newdata = np.array([(folder,q,numcharges,distance)+a+(lumoave,homoave,homolumoov)], dtype=data.dtype)
                data = append(data,newdata,axis = 0)
           
            

    
    
    data.sort()
    with open('./agglomeratedDFTresults.pic','wb') as anopenfile:
        pickle.dump(data,anopenfile)
        anopenfile.close()
    return data
    
def plotgroundstateDFTdata():
    """Plot results of workup ground  state DFT data from computations on Cd29S29 cluster with different numbers of point charges
    on it.  """   
    
    os.chdir('/home/chris/Dropbox/DataWeiss/160120/chargeresults')
    with open('./agglomeratedDFTresults.pic','rb') as anopenfile:
        data = pickle.load(anopenfile)
        anopenfile.close()
    fig=figure()
    ax = fig.add_subplot(131)
    plot_dfts(data,data['charge']/data['distance'],1240/data['gap'], {'numcharges':1024,'distance':0.5},'rs-')
    plot_dfts(data,data['charge']/data['distance'],1240/data['gap'], {'numcharges':1024,'distance':1},'ks-')
    plot_dfts(data,data['charge']/data['distance'],1240/data['gap'], {'numcharges':1024,'distance':2},'kx-')
    plot_dfts(data,data['charge']/data['distance'],1240/data['gap'], {'numcharges':1024,'distance':4},'kx-')
    plot([-16,-8,-2,-1,0,1,12,16],array([643.58748119,  535.6602877,   515.97869507,  513.79796138,  511.74115802,  509.93132377,  495.74221405,  492.20021435]),'rs', markersize = 10)
    
    plot_dfts(data,data['charge'],1240/data['gap'], {'charge':0},'ko',markersize = 5)
    ylabel('homo-lumo gap (nm)')
    xlabel('q/r')
    legend([ 'distance = 0.5 A','1 A','2 A','4 A','TDDFT 1A'])
    
    ax = fig.add_subplot(132)
    plot_dfts(data,data['charge']/data['distance'],data['homo'], {'numcharges':1024,'distance':1},'ro')
    plot_dfts(data,data['charge']/data['distance'],data['lumo'], {'numcharges':1024,'distance':1},'bo')
    plot_dfts(data,data['charge']/data['distance'],data['homo'], {'numcharges':1024,'distance':0.5},'rs')
    plot_dfts(data,data['charge']/data['distance'],data['lumo'], {'numcharges':1024,'distance':0.5},'bs')
    plot_dfts(data,data['charge']/data['distance'],data['homo'], {'charge':0},'ro',markersize = 10)
    plot_dfts(data,data['charge']/data['distance'],data['lumo'], {'charge':0},'bo',markersize = 10)

    ylabel('homo(lumo) energies (eV)')
    xlabel('q/r')
    legend(['distance = 1 A HOMO','distance = 1 A LUMO', 'distance = 0.5 A HOMO','distance = 0.5 A LUMO'])
    
    ax = fig.add_subplot(133)  ### number of charges
   
    plot_dfts(data,data['numcharges'],1240/data['gap'], {'charge':-1,'distance':0.5},'ks-')
    plot_dfts(data,data['numcharges'],1240/data['gap'], {'charge':-1,'distance':1},'rs-')
    
    xlabel('number of charges')
    ylabel('homo-lumo gap (nm)')
    annotate('total charge = -1\ndistance of charges from surface = 0.5 A', (0.1,0.8),xycoords = 'axes fraction')
    xlim(0,17)
    
    
    
    
    

    fig = figure()
    ax = fig.add_subplot(131)
    plot_dfts(data,data['distance'],data['energy'], {'numcharges':1024,'charge':-1},'rs')
    plot_dfts(data,data['distance'],data['energy'], {'numcharges':1024,'charge':-2},'bo')
    plot_dfts(data,data['distance'],data['energy'], {'numcharges':1024,'charge':-3},'kx')
    ylabel('DFT energies (eV)')
    xlabel('total charge (a. u.)')
    legend(['distance = 1 A HOMO','distance = 1 A LUMO', 'distance = 0.5 A HOMO','distance = 0.5 A LUMO'])
    
    ax = fig.add_subplot(132)  ### number of charges
    plot_dfts(data,data['charge']/data['distance'],data['lumoaveD'], {'numcharges':1024,'distance':4},'rs')
    plot_dfts(data,data['charge']/data['distance'],data['lumoaveD'], {'numcharges':1024,'distance':2},'rx')
    plot_dfts(data,data['charge']/data['distance'],data['lumoaveD'], {'numcharges':1024,'distance':1},'ro')
    plot_dfts(data,data['charge']/data['distance'],data['lumoaveD'], {'numcharges':1024,'distance':0.5},'r-')
    plot_dfts(data,data['charge']/data['distance'],data['homoaveD'], {'numcharges':1024,'distance':4},'bs')  
    plot_dfts(data,data['charge']/data['distance'],data['homoaveD'], {'numcharges':1024,'distance':2},'bx')
    plot_dfts(data,data['charge']/data['distance'],data['homoaveD'], {'numcharges':1024,'distance':1},'bo')
    plot_dfts(data,data['charge']/data['distance'],data['homoaveD'], {'numcharges':1024,'distance':0.5},'b-')
    ylabel('Average Radius of HOMO(LUMO) Density (Angstrom)')
    xlabel('q/r ')
    legend(['4A','2A', '1A','0.5 A'])
    
    ax = fig.add_subplot(133)  ### overlap of homo and lumo
    plot_dfts(data,data['charge']/data['distance'],data['homolumooverlap'], {'numcharges':1024,'distance':0.5},'s')
    plot_dfts(data,data['charge']/data['distance'],data['homolumooverlap'], {'numcharges':1024,'distance':1},'s')
    plot_dfts(data,data['charge']/data['distance'],data['homolumooverlap'], {'numcharges':1024,'distance':2},'s')
    plot_dfts(data,data['charge']/data['distance'],data['homolumooverlap'], {'numcharges':1024,'distance':4},'s')
    ylabel('Average Radius of HOMO(LUMO) Density (Angstrom)')
    xlabel('q/r')
    legend(['lumos 0.5 A','lumos 1A','lumos 2A','lumos 4A'])
    
    return 0
    
def plot_dfts(data,x,y,criteria,*args, **kwargs):
    """ special function for plotting the record arrays generated by the DFT workup"""
    
    crit = array([True]*data.size)
    
    for c in criteria:
        if c == 'zero':
            pass
        
        crit = np.logical_and(crit,data[c]==criteria[c])
    if 'zero' in criteria:
        crit[data['charge'==0]]=True
    order = np.argsort(x[crit])

    plot(x[crit][order],y[crit][order],*args, **kwargs)
    return None
    
def tddftresultsforchargedCd29S29():
    """process and plot the results of TDDFT data of Cd29S29 clusters with point charges"""
    
    fig1 = figure()
    ax1 = fig1.add_subplot(121)
    ax2 = fig1.add_subplot(122)
    figure()
    subplot(211)
    energies = array([])
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-16d1024r1_Cd29S29.log')))
    energies = np.append(energies, 1240/2.3149)#max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-8d1024r1_Cd29S29.log')))
    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-2d1024r1_Cd29S29.log')))
    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d1024r1_Cd29S29.log')))
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Cd29S29_LANL_B3LYP_TDDFT.log')))
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q1d1024r1_Cd29S29.log')))
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q2d1024r1_Cd29S29.log')))
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q3d1024r1_Cd29S29.log')))
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q12d1024r1_Cd29S29.log')))
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q16d1024r1_Cd29S29.log')))
    
    print energies
    
    ax1.plot([-16,-8,-2,-1,-0,1,2,3,12,16],energies,'s-')
    legend(['q=-16','-8','-2','-1','0','12'])
    
    subplot(212)
    energies = array([])
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d1024r1_Cd29S29.log')))
    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d16r1_Cd29S29.log')))
    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d8r1_Cd29S29.log')))
    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d4r1_Cd29S29.log')))
    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d2r1_Cd29S29.log')))
    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d1r1_Cd29S29.log')))
    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Cd29S29_LANL_B3LYP_TDDFT.log')))
    
    print energies
    ax2.plot([20,16,8,4,2,1],energies)
    legend(['n=1024','16','8','4','2','1'])
    return 0
    
def tddftabsorptionfigure():
    """process and plot the results of TDDFT data of Cd29S29 clusters with point charges"""
    
    fig1 = figure()
    ax1 = fig1.add_subplot(121)
    ax2 = fig1.add_subplot(122)

    energies = array([])
   # energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-16d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies = np.append(energies, min(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-8d1024r1_Cd29S29.log',lines=False,spectrum=True)))#1240/2.3149)#
    energies = np.append(energies,min(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-4d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies = np.append(energies,min(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-2d1024r1_Cd29S29.log',lines=False,spectrum=True)))
  #  energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies = np.append(energies, min(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Cd29S29_LANL_B3LYP_TDDFT.log',lines=False,spectrum=True)))
  #  energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q1d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies = np.append(energies, min(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q2d1024r1_Cd29S29.log',lines=False,spectrum=False)))
    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q3d1024r1_Cd29S29.log',lines=False,spectrum=False)))
    energies = np.append(energies, min(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q4d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies = np.append(energies, min(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q8d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q12d1024r1_Cd29S29.log',lines=False,spectrum=False)))
    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q16d1024r1_Cd29S29.log',lines=False,spectrum=False)))
    
    
    energies2 = array([])
   # energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-16d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies2 = np.append(energies2, NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-8d1024r1_Cd29S29.log',lines=False,spectrum=False)[1])#1240/2.3149)#
    energies2 = np.append(energies2,NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-4d1024r1_Cd29S29.log',lines=False,spectrum=False)[1])
    energies2 = np.append(energies2,NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-2d1024r1_Cd29S29.log',lines=False,spectrum=False)[1])
  #  energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies2 = np.append(energies2, NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Cd29S29_LANL_B3LYP_TDDFT.log',lines=False,spectrum=False)[1])
  #  energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q1d1024r1_Cd29S29.log',lines=False,spectrum=True)))
    energies2 = np.append(energies2, NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q2d1024r1_Cd29S29.log',lines=False,spectrum=False)[1])
    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q3d1024r1_Cd29S29.log',lines=False,spectrum=False)))
    energies2 = np.append(energies2, NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q4d1024r1_Cd29S29.log',lines=False,spectrum=False)[1])
    energies2 = np.append(energies2, NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q8d1024r1_Cd29S29.log',lines=False,spectrum=False)[1])
    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q12d1024r1_Cd29S29.log',lines=False,spectrum=False)))
    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Q16d1024r1_Cd29S29.log',lines=False,spectrum=False)))
    
    
    ax2.set_ylabel('absorbance')
    ax2.set_xlabel('energy (eV)')
    
    ax1.plot([-8,-4,-2,0,2,4,8],energies,'s-')
    ax1.plot([-8,-4,-2,0,2,4,8],energies2,'s-')
    legend(list(str(i) for i in [-8,-4,-2,0,2,4,8]))
    ax1.set_ylabel(' exciton energy (eV)')
    ax1.legend(['first exciton', 'second exciton'])
    
    ax1.set_xlabel('magnitude of charge on the surface (a. u.)')
    print energies2-energies
 
#    figure()
#    ax1 = fig1.add_subplot(121)
#    ax2 = fig1.add_subplot(122)
#    
#    energies = array([])
#    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d1024r1_Cd29S29.log')))
#    energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d16r1_Cd29S29.log')))
#    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d8r1_Cd29S29.log')))
#    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d4r1_Cd29S29.log')))
#    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d2r1_Cd29S29.log')))
#    energies = np.append(energies,max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/tddftQ-1d1r1_Cd29S29.log')))
#    #energies = np.append(energies, max(NWchemTDDFT('/home/chris/Dropbox/DataWeiss/160120/tddftcharges_results/Cd29S29_LANL_B3LYP_TDDFT.log')))
#    
#    print energies
#    ax2.plot([20,16,8,4,2,1],energies)
#    legend(['n=1024','16','8','4','2','1'])
    return 0

def averagefieldatsurface():
    """calculate the average energy of a point charge at the surface of a sphere which is surrounded by many point charges"""
    
    ## a point on the surface at [0,0,0.7] feels the effect of 1024 charges around the
  ### concentric sphere, with diameter r  
    average1overRs = array([])
    d = 0  ## radius of nanoparticle
    xs = arange(1,50,0.1)
    totalchargeofsurface=1.0
    numberofpointcharges = 2
    for r in xs :  ### radius of the charge sphere
        z = pointsonsphere(numberofpointcharges,r)
        if any(np.sqrt(sum(z[:]**2,axis = 1))-r>0.1):
            print 'error'
            return -1
        
        average1overR =(totalchargeofsurface/numberofpointcharges)*1E10*(z[:]-array([0,0,d]))/transpose([sum((z[:]-array([0,0,d]))**2,axis = 1)])**(3/2)   ### units m^-1
        average1overR = sum(average1overR,axis = 0)
        average1overRs = append(average1overRs, sqrt(sum(average1overR**2)))
    eps = 9
    eps0 = 8.854E-12
    inkilojoulespermol = 6.02E23*average1overRs*1.602E-19**2/(4*pi*eps*eps0)/1000
    plot(xs, inkilojoulespermol)
    xlabel("distance of point charges from surface Angstroms")
    ylabel("total coulomb energy (kJ/mol)")
    return 0