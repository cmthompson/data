# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 22:57:17 2015

@author: chris
"""
#from __future__ import division
from numpy import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import itertools
import numpy


def makeacluster(radius,charges = 0, distance =4,cadmium_term=False,
                 module = 'nwchem',
                 removehangeroners = True, surfacepassivation=True,_plot=True):
    """ Prepare an xyzfile of cluster of CdS or CdSe with a specified radius and a certain number of point charges around it"""

    
    #########cell vectors in angstroms
    a = array([4.16,0,0])
    b = 4.16*array([cos(radians(60)),sin(radians(60)),0])
    c = array([0,0,6.756])
    
    
    ###### cell vectors in angstroms for CdSe
    a = array([4.2985,0,0])
    b = 4.2985*array([cos(radians(60)),sin(radians(60)),0])
    c = array([0,0,7.0152])
        
        
    cellvectors = array([a,b,c])
    
    
    
    ##Atomic positions in fractional coordinates in cell
    ### For Wurtzite#####
    Cdcoord = array([[0,0,0],[0.33,0.333,0.5]])
    Scoord = array([[0,0,0.375],[0.3333,0.3333,0.875]])
    tetarray = array([[0,0,0.375],[-0.333,-0.333,-0.125],[0.6666,-0.3333,-0.125],[-0.3333,0.6666,-0.125]])
    Cdcoord= dot(Cdcoord,cellvectors)
    Scoord=dot(Scoord,cellvectors)
  

    Cdcoordlist1 = ndarray((0,3))
    Cdcoordlist2 = ndarray((0,3))
    Scoordlist1 = ndarray((0,3))
    Scoordlist2 = ndarray((0,3))
    OHcoordlist=ndarray((0,3))
    Hcoordlist=ndarray((0,3))
    
    numbertorepeat = int(radius/b[1])+2
    
    for p in itertools.product(range(-numbertorepeat,numbertorepeat+1),repeat=3):
        
        x= dot(array([p]),cellvectors)+Cdcoord[0]
        x2=dot(array([p]),cellvectors)+Cdcoord[1]
        y = dot(array([p]),cellvectors)+Scoord[0]
        y2 = dot(array([p]),cellvectors)+Scoord[1]
        
        Cdcoordlist1 = append(Cdcoordlist1,x,axis=0)
        Cdcoordlist2 = append(Cdcoordlist2,x2,axis=0)
        Scoordlist1 = append(Scoordlist1,y,axis=0)
        Scoordlist2 = append(Scoordlist2,y2,axis=0)
    
    
    
    tetcoord=array([[ -2.08208000e+00,   1.19968767e+00,  -8.44500000e-01],
                     [ -2.08000000e-03,  -2.40297801e+00 , -8.44500000e-01],
                     [  2.07792000e+00,   1.19968767e+00,  -8.44500000e-01],
                     [  0.00000000e+00,   0.00000000e+00,   2.53350000e+00]])#ndarray((0,3))

    tetcoord2 = array( [[ -2.06544000e+00,  -1.19968767e+00,  -8.44500000e-01],
                             [  1.43520000e-02,   1.08079970e-03,   2.53350000e+00],
                             [  1.45600000e-02,   2.40297801e+00,  -8.44500000e-01],
                             [  2.09456000e+00,  -1.19968767e+00,  -8.44500000e-01]])
    Cdtetcoord1=tetcoord*array([[1,1,-1]])
    Cdtetcoord2=tetcoord2*array([[1,1,-1]])
    
    extrashift = (a+b+c)/2

    Cdcoordlist1[:]+=extrashift#array([[radius+extrashift,radius+extrashift,radius+extrashift]])
    Scoordlist1[:]+=extrashift#array([[radius+extrashift,radius+extrashift,radius+extrashift]])
    Cdcoordlist2[:]+=extrashift#array([[radius+extrashift,radius+extrashift,radius+extrashift]])
    Scoordlist2[:]+=extrashift#array([[radius+extrashift,radius+extrashift,radius+extrashift]])
    
    x= where(numpy.sum(Cdcoordlist1**2,axis=1)<radius**2)[0]
    Cdcoordlist1=Cdcoordlist1[x]
    
    x= where(numpy.sum(Cdcoordlist2**2,axis=1)<radius**2)[0]
    Cdcoordlist2=Cdcoordlist2[x]
   
    
    
    Scoordlist1=Scoordlist1[numpy.sum(Scoordlist1**2,axis=1)<radius**2]
    Scoordlist2=Scoordlist2[numpy.sum(Scoordlist2**2,axis=1)<radius**2]
    
   
    extrashift = 10
    Cdcoordlist1[:]+=array([[radius+extrashift,radius+extrashift,radius+extrashift]])
    Scoordlist1[:]+=array([[radius+extrashift,radius+extrashift,radius+extrashift]])
    Cdcoordlist2[:]+=array([[radius+extrashift,radius+extrashift,radius+extrashift]])
    Scoordlist2[:]+=array([[radius+extrashift,radius+extrashift,radius+extrashift]])
                  

    print '--------------'


    if cadmium_term:
        """Cadmium Termination.  Still not working"""
        ### TODO implement this.  
        print 'terminating the surface with cadmium atoms'

                       
        for i in Scoordlist1:
            remove = array([])
            bondlength=1.83  ## angstroms
            x= append(Cdcoordlist1,Cdcoordlist2,axis=0)[numpy.sum((append(Cdcoordlist1,Cdcoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==4:
                pass
            else:
                
                for s in Cdtetcoord1:
                    if not any(sum((x[:]-s)**2,axis=1)<0.5):
                       print "    adding Cd1 Atom at",[i+s]
                       Cdcoordlist1= append(Cdcoordlist1,[i+s],axis=0)
             
        for i in Scoordlist2:
            bondlength=1.83 ## angstroms
            x= append(Cdcoordlist1,Cdcoordlist2,axis=0)[numpy.sum((append(Cdcoordlist1,Cdcoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==4:
                pass
            else:
                #print 'undercoordinated atom found'
                for s in Cdtetcoord2:
                    if not any(sum((x[:]-s)**2,axis=1)<0.5):
                       print "    adding Cd2 Atom at",[i+s]
                       Cdcoordlist2= append(Cdcoordlist2,[i+s],axis=0)
    
        
                       
    if removehangeroners:
        """removing Cd or S/Se atoms connected by only one bond"""
        print '--------------'  
        print "removing core atoms connected by only one bond"
        removefromCd1 = []
        removefromCd2 = []
        removefromS1 = []
        removefromS2 = []
        removefromS2_2 = array([])
        removefromCd1_2 = array([])
        for i in Cdcoordlist1:
        
            x= append(Scoordlist1,Scoordlist2,axis=0)[numpy.sum((append(Scoordlist1,Scoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
    
            if x.shape[0]==1:
                        removefromCd1 = append(removefromCd1,where(all(Cdcoordlist1==i,axis=1))[0])
        

        for i in Cdcoordlist2:
            x= append(Scoordlist1,Scoordlist2,axis=0)[numpy.sum((append(Scoordlist1,Scoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==1:
                removefromCd2 = append(removefromCd2,where(all(Cdcoordlist2==i,axis=1))[0])
                
                          
                           
            
        for i in Scoordlist1:
            x= append(Cdcoordlist1,Cdcoordlist2,axis=0)[numpy.sum((append(Cdcoordlist1,Cdcoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==1:
                removefromS1 = append(removefromS1,where(all(Scoordlist1==i,axis=1))[0])
            
        
               
        for i in Scoordlist2:
            x= append(Cdcoordlist1,Cdcoordlist2,axis=0)[numpy.sum((append(Cdcoordlist1,Cdcoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==1:
                removefromS2 = append(removefromS2,where(all(Scoordlist2==i,axis=1))[0])
                
        print "    removed", len(removefromCd1), "atoms from Cd1"
        print "    removed", len(removefromCd2), "atoms from Cd2"
        print "    removed", len(removefromS1), "atoms from S1"
        print "    removed", len(removefromS2), "atoms from S2"
        Cdcoordlist1 = delete(Cdcoordlist1,removefromCd1,axis=0)  
        Cdcoordlist2 = delete(Cdcoordlist2,removefromCd2,axis=0)  
        Scoordlist1 = delete(Scoordlist1,removefromS1,axis=0)  
        Scoordlist2 = delete(Scoordlist2,removefromS2,axis=0)  
        

    listofcoordlists = [Cdcoordlist1,Cdcoordlist2,Scoordlist1,Scoordlist2,OHcoordlist,Hcoordlist]
    alllists = ndarray((0,3))                                                                   


    """Surface passivation with cations and anions"""
    if surfacepassivation:
        CdOHbondlength =2.7
        for i in Cdcoordlist1:
            bondlength=CdOHbondlength
            x= append(Scoordlist1,Scoordlist2,axis=0)[numpy.sum((append(Scoordlist1,Scoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            
            if x.shape[0]==4:
                pass
            else:
                #print 'undercoordinated atom found'
                for s in tetcoord:
                    if not any(sum((x[:]-s)**2,axis=1)<0.5):# and not any(sum((OHcoordlist[:]-(i+s/sqrt(sum(s**2))*bondlength))**2,axis=1)<0.5):
                       
                        OHcoordlist= append(OHcoordlist,[i+s/sqrt(sum(s**2))*bondlength],axis=0)
                    else:
                         pass
                        # print min(sum((OHcoordlist[:]-(i+s/sqrt(sum(s**2))*bondlength))**2,axis=1))
    
        for i in Cdcoordlist2:
            bondlength=CdOHbondlength ## angstroms
            x= append(Scoordlist1,Scoordlist2,axis=0)[numpy.sum((append(Scoordlist1,Scoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==4:
                pass
            else:
                #print 'undercoordinated atom found'
                for s in tetcoord2:
                    if not any(sum((x[:]-s)**2,axis=1)<0.5):# and not any(sum((OHcoordlist[:]-(i+s/sqrt(sum(s**2))*bondlength))**2,axis=1)<0.5):
                        
                        OHcoordlist= append(OHcoordlist,[i+s/sqrt(sum(s**2))*bondlength],axis=0)
                    else:
                        pass#print min(sum((OHcoordlist[:]-(i+s/sqrt(sum(s**2))*bondlength))**2,axis=1))
    
        for i in Scoordlist1:
            remove = array([])
            bondlength=1.83  ## angstroms
            x= append(Cdcoordlist1,Cdcoordlist2,axis=0)[numpy.sum((append(Cdcoordlist1,Cdcoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==4:
                pass
            else:
                #print 'undercoordinated atom found'
                for s in Cdtetcoord1:
                    if not any(sum((x[:]-s)**2,axis=1)<0.5) and not any(sum((Hcoordlist[:]-(i+s/sqrt(sum(s**2))*bondlength))**2,axis=1)<2):
                       
                       Hcoordlist= append(Hcoordlist,[i+s/sqrt(sum(s**2))*bondlength],axis=0)
             
        for i in Scoordlist2:
            bondlength=1.83 ## angstroms
            x= append(Cdcoordlist1,Cdcoordlist2,axis=0)[numpy.sum((append(Cdcoordlist1,Cdcoordlist2,axis=0)[:]-i)**2,axis=1)<9]-i
            if x.shape[0]==4:
                pass
            else:
                #print 'undercoordinated atom found'
                for s in Cdtetcoord2:
                    if not any(sum((x[:]-s)**2,axis=1)<0.5) and not any(sum((Hcoordlist[:]-(i+s/sqrt(sum(s**2))*bondlength))**2,axis=1)<2):
                       
                       Hcoordlist= append(Hcoordlist,[i+s/sqrt(sum(s**2))*bondlength],axis=0)
    
        
    
        listofcoordlists = [Cdcoordlist1,Cdcoordlist2,Scoordlist1,Scoordlist2,OHcoordlist,Hcoordlist]
     
     
     
    ###########Remove Duplicate Atoms 
    for i in listofcoordlists[0:4]:
       
        OHcoordlist = delete(OHcoordlist,checkduplicates(i,OHcoordlist)[1],axis=0)
        Hcoordlist = delete(Hcoordlist,checkduplicates(i,Hcoordlist)[1],axis=0)
    listofcoordlists = [Cdcoordlist1,Cdcoordlist2,Scoordlist1,Scoordlist2,OHcoordlist,Hcoordlist]
    
    for i in range(len(listofcoordlists)):
        for j in range(i,len(listofcoordlists)):
            
            
            x = checkduplicates(listofcoordlists[i],listofcoordlists[j],threshold=1)
           
            if x[0].size>0:
                if j==i:
#                    print i,j,x
#                    print listofcoordlists[i][x[0]]
#                    print  listofcoordlists[i][x[1]]
                    listofcoordlists[i]=delete(listofcoordlists[i],x[1],axis=0)
                    
                    
                    print 'duplicate atoms found in same list! THE DUPLICATE ATOMS WERE REMOVED'
                else:
                    print i,j,x
                    print 'duplicate atoms found! ERROR YOU NEED TO FIX THIS'
   
          
    (Cdcoordlist1,Cdcoordlist2,Scoordlist1,Scoordlist2,OHcoordlist,Hcoordlist)=listofcoordlists[:]

 
  ####################3 add charges aroundthe particls
    chargecoordlist=pointsonsphere(charges,radius+distance)+array([[radius,radius,radius]])

   ###############################    
    if _plot: 
        fig=figure()
        ax = fig.add_subplot(111, projection='3d')
            
                
        Axes3D.scatter(ax,Cdcoordlist1[:,0],Cdcoordlist1[:,1],Cdcoordlist1[:,2],c='b')
        Axes3D.scatter(ax,Cdcoordlist2[:,0],Cdcoordlist2[:,1],Cdcoordlist2[:,2],c='b')
        Axes3D.scatter(ax,Scoordlist1[:,0],Scoordlist1[:,1],Scoordlist1[:,2],c='r') 
        Axes3D.scatter(ax,Scoordlist2[:,0],Scoordlist2[:,1],Scoordlist2[:,2],c='r')    
        Axes3D.scatter(ax,OHcoordlist[:,0],OHcoordlist[:,1],OHcoordlist[:,2],c='y')
        Axes3D.scatter(ax,Hcoordlist[:,0],Hcoordlist[:,1],Hcoordlist[:,2],c='k')
        Axes3D.scatter(ax,chargecoordlist[:,0],chargecoordlist[:,1],chargecoordlist[:,2],c='k')



    numcoreatoms = Cdcoordlist1.shape[0]*1+Cdcoordlist2.shape[0]+Scoordlist1.shape[0]*1+Scoordlist2.shape[0]*1
    print 'cadmium atoms:',Cdcoordlist1.shape[0]*1+Cdcoordlist2.shape[0]
    print 'sulfur atoms:', Scoordlist1.shape[0]*1+Scoordlist2.shape[0]
    print 'negative cations:', OHcoordlist.shape[0]
    print 'positive cations:', Hcoordlist.shape[0]
    print 'number of core atoms', numcoreatoms
    numatoms= Cdcoordlist1.shape[0]*1+Cdcoordlist2.shape[0]+Scoordlist1.shape[0]*1+Scoordlist2.shape[0]*1+OHcoordlist.shape[0]*1+Hcoordlist.shape[0]*1
    
    print 'num electrons with ECP:', str(10*(Cdcoordlist1.shape[0]+Cdcoordlist2.shape[0])+
                                                8*(Scoordlist1.shape[0]+Scoordlist2.shape[0])+
                                                OHcoordlist.shape[0]*10+
                                                Hcoordlist.shape[0]*0)

    print 'num atoms:', numatoms, 'total charge:', str(Cdcoordlist1.shape[0]*2+Cdcoordlist2.shape[0]*2+
                                                Scoordlist1.shape[0]*-2+Scoordlist2.shape[0]*-2+
                                                OHcoordlist.shape[0]*-1+
                                                Hcoordlist.shape[0]*1)
    
    commentline = str('cadmium atoms: '+str(Cdcoordlist1.shape[0]*1+Cdcoordlist2.shape[0])+
                ', sulfur atoms: '+ str(Scoordlist1.shape[0]*1+Scoordlist2.shape[0])+
                ', total charge of cluster: '+ str(Cdcoordlist1.shape[0]*2+
                  +Cdcoordlist2.shape[0]*2
                                +Scoordlist1.shape[0]*-2
                                +Scoordlist2.shape[0]*-2
                                +OHcoordlist.shape[0]*-1
                                +Hcoordlist.shape[0]*1)+
                ', number of point charges:'+str(charges)+
                ' removehangeroners= '+str(removehangeroners)
                +'\n')
    with open('/home/chris/Desktop/CdSParticle.xyz','wb') as f:#'+str(numcoreatoms)+'core'+str(charges)+'charges.xyz','wb') as f:
        f.write(str(numatoms)+'\n')
        f.write(commentline)
        
        
        for i in append(Cdcoordlist1,Cdcoordlist2,axis=0):
            f.write('Cd '+str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
        for i in append(Scoordlist1,Scoordlist2,axis=0):
            f.write('S '+str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
        for i in OHcoordlist:
            f.write('Cl '+str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
        for i in Hcoordlist:
            f.write('H '+str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
        for i in chargecoordlist:
            
            if module == 'nwchem':
                f.write('bq1 '+str(i[0])+' '+str(i[1])+' '+str(i[2])+' charge -1\n')
            elif module == 'orca':
                f.write('Q -1'+str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
        f.close()
    return None


 
def pointsonsphere(n,rad):
    if n ==0:
        return ndarray((0,3))
    
 
    golden_angle = numpy.pi * (3 - numpy.sqrt(5))
    theta = golden_angle * numpy.arange(n)
    z = numpy.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
    radius = numpy.sqrt(1 - z * z)
     
    points = numpy.zeros((n, 3))
    points[:,0] = radius * numpy.cos(theta)
    points[:,1] = radius * numpy.sin(theta)
    points[:,2] = z
    
    points = numpy.zeros((n, 3))
    
    points[:,0] = radius * numpy.cos(theta)*rad
    points[:,1] = radius * numpy.sin(theta)*rad
    points[:,2] = z*rad
    
   # ax = fig.add_subplot(111, projection='3d')
    
    return points
    
def perm(n, seq):
    coordlist = ndarray((0,3))
    for p in itertools.product(seq, repeat=n):
        coordlist = append(coordlist,array([p]),axis=0)
    return coordlist
    
def rotate(vector, aroundvector, bydegrees):
    if all(vector == 0):
        return vector
    bydegrees = radians(bydegrees)
    apar = dot(vector,aroundvector)*aroundvector/sum(aroundvector**2)**2
    
    aper = vector - apar
    
    w = numpy.cross(aroundvector,aper)
    
    aperrot = sqrt(sum(aper**2))*(cos(bydegrees)*aper/sqrt(sum(aper**2))+sin(bydegrees)*w/sqrt(sum(w**2)))
    return aperrot+apar

def symmetrizeC3(geometry,angle):
     ##### Now symmetrize to C3V by taking the average of the rotated clusters
    center = mean(geometry,axis=0)
    geometry[:]-=center
    
    geometryrot = ndarray(geometry.shape)
    geometryrot2 = ndarray(geometry.shape)
   
    for i in geometry:
        rotatedatom = rotate(i,array([0,0,1]),120)
        if min(sum((geometry[:]-array([rotatedatom]))**2, axis = 1))>0.01:
            print "symmetry error detected"
      
        closestatom = argmin(sum((geometry[:]-array([rotatedatom]))**2, axis = 1)) 
               
        geometryrot[closestatom] = rotatedatom
    
    for i in geometry:
        rotatedatom = rotate(i,array([0,0,1]),240)
        
        closestatom = argmin(sum((geometry[:]-array([rotatedatom]))**2, axis = 1)) 
               
        geometryrot2[closestatom] = rotatedatom
        
    geometry_mean = (geometry+geometryrot+geometryrot2)/3    
    
    
   
    with open('/home/chris/Desktop/rotatedgeo.xyz', 'wb') as f:
        f.write('66\nrotetateclsuter\n')
        for i in range(66):
            if i<33:
                f.write('Cd '+str(geometry_mean[i][0])+' '+str(geometry_mean[i][1])+' '+str(geometry_mean[i][2])+'\n')
            else:
                f.write('S '+str(geometry_mean[i][0])+' '+str(geometry_mean[i][1])+' '+str(geometry_mean[i][2])+'\n')
        f.close()
    
    print geometry_mean-geometry
    print checkduplicates(geometry_mean, geometry_mean)
    return 0
        

def addMePA(xyzfile):
    allatoms = array([[ 3.73479929, 4.70148859, 1.08453093],
    [ -0.45777444, 5.46389693, 1.4292281],
    [ 3.99256744, 1.42776590, 0.17101104],
    [ -3.23253881, 2.73557278, 0.17389669],
    [ 0.00277074, -0.00738694, -0.25202305],
    [ 0.00013075, -0.01651818, 6.38550475],
    [ 4.95224596, -2.34187282, 1.40722903],
    [ -5.94120956, 0.87524543, 1.07437676],
    [ -4.49880283, -3.12964705, 1.41125120],
    [ -0.75991343, -4.18114040, 0.15674939],
    [ 2.20650584, -5.59236598, 1.05533723],
    [ 0.96252338, 4.36198028, -1.42717920],
    [ 4.81764990, 3.19485539, -2.77622069],
    [ -1.73138701, 4.16471485, -3.46626689],
    [ -2.55306625, 2.91594904, 4.51767596],
    [ 1.50206160, 1.89218763, -4.09811156],
    [ 0.98329346, 2.34564526, 3.59745292],
    [ 4.46454220, -0.56369361, -3.4747580],
    [ 3.80849139, 0.73094006, 4.51682402],
    [ -5.17525966, 2.60395011, -2.77772471],
    [ -2.37975873, 0.36837952, -4.07698444],
    [ -2.53045715, -0.33540920, 3.59449188],
    [ 0.88144546, -2.23348348, -4.09931249],
    [ 1.54175224, -2.03600240, 3.58132334],
    [ 3.31038688, -3.02228474, -1.45695989],
    [ -4.27528734, -1.35796536, -1.44852313],
    [ -2.72513737, -3.56991001, -3.48209886],
    [ -1.25794100, -3.68787964, 4.50492503],
    [ 0.35261163, -5.75579556, -2.79552608],
    [ 0.74123700, 4.43697518, -4.02956973],
    [ -0.64622546, 4.44077927, 3.74595547],
    [ 4.18920378, 1.70873805, -4.66166552],
    [ 3.47773976, 2.51701231, 2.54165082],
    [ -3.55947851, 2.79106662, -4.65228042],
    [ -3.92545305, 1.73860005, 2.54076867],
    [ -0.00271263, 0.01013464, -5.02518153],
    [ -0.00332939, -0.00785361, 2.48730355],
    [ 3.46776368, -2.83946171, -4.05736012],
    [ 4.17115409, -1.67896151, 3.73089394],
    [ -4.20196677, -1.57060486, -4.04933515],
    [ -3.53060015, -2.79614465, 3.73504438],
    [ -0.62463325, -4.46372023, -4.67635054],
    [ 0.44748177, -4.28341329, 2.52481152],
    [ 1.79334057, 6.15314947, 0.33208226],
    [ 5.55303831, 3.65754425, -0.46937522],
    [ -1.77455171, 4.68869481, -0.74204629],
    [ 1.74568849, 1.79694051, -1.26045912],
    [ 1.94632660, 1.65173243, 6.00004182],
    [ 4.95599592, -0.80678637, -0.75922862],
    [ -5.94441654, 2.98322195, -0.46759778],
    [ -2.43853084, 0.59877699, -1.25432294],
    [ -2.41524712, 0.83987083, 5.99734552],
    [ 0.69721678 ,-2.42073954, -1.26854049],
    [ 0.46924459, -2.53096433, 5.98775239],
    [ 4.43666586, -4.63906404, 0.30520008],
    [ -6.23050849, -1.53090073, 0.31876024],
    [ -3.17891818, -3.89146035, -0.76368896],
    [ 0.38866467, -6.64405085, -0.49626228]])
    atomcenter = sum(allatoms,axis=0)
    atomcenter = array([0,0,0])
    
#    C       -1.4255400000      2.275800,0000      0.1498900000                         
#    P       -2.0727700000      0.8165500000     -0.7496400000 
#    H       -0.3427600000      2.2060900000      0.0812300000                       
#    H       -1.8434700000      3.0857900000     -0.4746900000                 
#    H       -1.8709200000      2.2658600000      1.1464000000                 
#    O       -1.4521900000      1.1860200000     -2.1113400000                 
#    O       -3.5385700000      1.1054100000     -0.4997000000                 
#    O       -1.5251800000     -0.4706000000     -0.1678200000  

    MePAvecs = array([[  -1.4255400000,2.2758000000,0.1498900000],                  
                [ -2.0727700000 ,0.8165500000,-0.7496400000 ],
                [ -0.3427600000 ,     2.2060900000,      0.0812300000  ],                     
                [ -1.8434700000 ,     3.0857900000,     -0.4746900000   ],              
                 [-1.8709200000 ,     2.2658600000,      1.1464000000  ],               
                [ -1.4521900000 ,     1.1860200000,     -2.1113400000   ],              
                 [-3.5385700000 ,     1.1054100000,     -0.4997000000  ],               
                  [ -1.5251800000,     -0.4706000000,     -0.1678200000 ]])
 
    MePAvecs-=MePAvecs[1]  ### put phosphorous at center
    CPvec = MePAvecs[0]
    Cdtarget = array([-1.2579410000  ,   -3.6878796400  ,    4.5049250300 ])  
    Pposition = (Cdtarget)/sqrt(sum((Cdtarget)**2)) 
    
     
    rotationaround = numpy.cross(CPvec,Pposition)
    #bydegrees = np.degrees(np.arcsin(sqrt(sum(rotationaround**2))/sqrt(sum(CPvec**2))/sqrt(sum(Pposition**2))))
    bydegrees=np.degrees(np.arccos(np.dot(CPvec,Pposition)/sqrt(np.dot(CPvec,CPvec)*np.dot(Pposition,Pposition))))
    print bydegrees    
    rotatedvecs = np.ndarray((0,3))    

    print '___________________'
    for i in MePAvecs:
        print i
        rotatedvecs = append(rotatedvecs, [rotate(i,rotationaround,bydegrees)+3*Pposition + Cdtarget],axis = 0) 
    print rotatedvecs
    #rotatedvecs[:]+=3*Pposition + Cdtarget
    print numpy.cross(rotatedvecs[0],Pposition)
     
#    return  0
#    finalP  = 3*Pposition + Cdtarget + atomcenter
#    v = numpy.cross(CPvec, Pposition)
#    s = sqrt(sum(v**2))#/sqrt(sum(CPvec**2))/sqrt(sum(Pposition**2))
#
#    c =   numpy.dot(CPvec, Pposition)#/sqrt(sum(CPvec**2))/sqrt(sum(Pposition**2))
#    print degrees(arcsin(s)), degrees(arccos(c))
#    if c<0:
#        c+=pi
#    vx = array([[0,-v[2],v[1]],
#                [v[2],0,-v[0]],
#                [-v[1],v[0],0]])
#    I = array([[1,0,0],[0,1,0],[0,0,1]])
#    rotationmatrix = I + vx + np.dot(vx,vx)*(1-c)/s**2
#  
#    rotatedvecs =transpose(np.dot(rotationmatrix,transpose(MePAvecs)))
#    #shift= rotatedvecs[1]-wherePshouldendup
#    rotatedvecs[:]+=finalP
#   # rotatedvecs[:]+=atomcenter#-=shift
    with open('/home/chris/Desktop/addMePa.xyz','wb') as f:#'+str(numcoreatoms)+'core'+str(charges)+'charges.xyz','wb') as f:
        f.write('66\n')
        f.write('MePA\n')
        f.write('''Cd 3.73479929 4.70148859 1.08453093
Cd -0.45777444 5.46389693 1.42922813
Cd 3.99256744 1.42776590 0.17101104
Cd -3.23253881 2.73557278 0.17389669
Cd 0.00277074 -0.00738694 -0.25202305
Cd 0.00013075 -0.01651818 6.38550475
Cd 4.95224596 -2.34187282 1.40722903
Cd -5.94120956 0.87524543 1.07437676
Cd -4.49880283 -3.12964705 1.41125120
Cd -0.75991343 -4.18114040 0.15674939
Cd 2.20650584 -5.59236598 1.05533723
Cd 0.96252338 4.36198028 -1.42717920
Cd 4.81764990 3.19485539 -2.77622069
Cd -1.73138701 4.16471485 -3.46626689
Cd -2.55306625 2.91594904 4.51767596
Cd 1.50206160 1.89218763 -4.09811156
Cd 0.98329346 2.34564526 3.59745292
Cd 4.46454220 -0.56369361 -3.47475802
Cd 3.80849139 0.73094006 4.51682402
Cd -5.17525966 2.60395011 -2.77772471
Cd -2.37975873 0.36837952 -4.07698444
Cd -2.53045715 -0.33540920 3.59449188
Cd 0.88144546 -2.23348348 -4.09931249
Cd 1.54175224 -2.03600240 3.58132334
Cd 3.31038688 -3.02228474 -1.45695989
Cd -4.27528734 -1.35796536 -1.44852313
Cd -2.72513737 -3.56991001 -3.48209886
Br -1.25794100 -3.68787964 4.50492503
Cd 0.35261163 -5.75579556 -2.79552608
S 0.74123700 4.43697518 -4.02956973
S -0.64622546 4.44077927 3.74595547
S 4.18920378 1.70873805 -4.66166552
S 3.47773976 2.51701231 2.54165082
S -3.55947851 2.79106662 -4.65228042
S -3.92545305 1.73860005 2.54076867
S -0.00271263 0.01013464 -5.02518153
S -0.00332939 -0.00785361 2.48730355
S 3.46776368 -2.83946171 -4.05736012
S 4.17115409 -1.67896151 3.73089394
S -4.20196677 -1.57060486 -4.04933515
S -3.53060015 -2.79614465 3.73504438
S -0.62463325 -4.46372023 -4.67635054
S 0.44748177 -4.28341329 2.52481152
S 1.79334057 6.15314947 0.33208226
S 5.55303831 3.65754425 -0.46937522
S -1.77455171 4.68869481 -0.74204629
S 1.74568849 1.79694051 -1.26045912
S 1.94632660 1.65173243 6.00004182
S 4.95599592 -0.80678637 -0.75922862
S -5.94441654 2.98322195 -0.46759778
S -2.43853084 0.59877699 -1.25432294
S -2.41524712 0.83987083 5.99734552
S 0.69721678 -2.42073954 -1.26854049
S 0.46924459 -2.53096433 5.98775239
S 4.43666586 -4.63906404 0.30520008
S -6.23050849 -1.53090073 0.31876024
S -3.17891818 -3.89146035 -0.76368896
S 0.38866467 -6.64405085 -0.49626228\n''')
        f.write('C '+str(rotatedvecs[0,0])+' '+str(rotatedvecs[0,1])+' '+str(rotatedvecs[0,2])+'\n')
        f.write('P '+str(rotatedvecs[1,0])+' '+str(rotatedvecs[1,1])+' '+str(rotatedvecs[1,2])+'\n')
        f.write('H '+str(rotatedvecs[2,0])+' '+str(rotatedvecs[2,1])+' '+str(rotatedvecs[2,2])+'\n')
        f.write('H '+str(rotatedvecs[3,0])+' '+str(rotatedvecs[3,1])+' '+str(rotatedvecs[3,2])+'\n')
        f.write('H '+str(rotatedvecs[4,0])+' '+str(rotatedvecs[4,1])+' '+str(rotatedvecs[4,2])+'\n')
        f.write('O '+str(rotatedvecs[5,0])+' '+str(rotatedvecs[5,1])+' '+str(rotatedvecs[5,2])+'\n')
        f.write('O '+str(rotatedvecs[6,0])+' '+str(rotatedvecs[6,1])+' '+str(rotatedvecs[6,2])+'\n')
        f.write('O '+str(rotatedvecs[7,0])+' '+str(rotatedvecs[7,1])+' '+str(rotatedvecs[7,2])+'\n')
        f.close()
    return rotatedvecs
def checkduplicates(a,b,threshold=0.001):
    infirstarray = array([],dtype = int)
    insecondarray = array([],dtype=int)
    if all(a==b):
        for x in range(a.shape[0]):

             r= where(sum((b[:][x+1:]-a[x])**2, axis = -1)<threshold)[0]+x+1
             if r.size>2:
                 print "ERROR IN CHECK DUPLICATES.  TRIPLICATE FOUND"
             
             insecondarray=numpy.append(insecondarray,r)
             if r.size>0:
                 
                 infirstarray=numpy.append(infirstarray,x)
    
    
    else:
        for x in range(a.shape[0]):
             #print x
             
             r= where(sum((b[:]-a[x])**2, axis = -1)<threshold)[0]
             if r.size>2:
                 print "ERROR IN CHECK DUPLICATES.  TRIPLICATE FOUND"
             #print r
             insecondarray=numpy.append(insecondarray,r)
             if r.size>0:
                 
                 infirstarray=numpy.append(infirstarray,x)
    return (infirstarray,insecondarray)
    


    
def prepdftcharged(charges,chargemag,totalcharge,d=0):
    """used to create input files for nwchem for Cd29S29 clusters with point charges"""
    z=array([[3.73479929,4.70148859,1.08453093],
    [-0.45777444,5.46389693,1.42922813],
    [3.99256744,1.42776590,0.17101104],
    [-3.23253881,2.73557278,0.17389669],
    [0.00277074,-0.00738694,-0.25202305],
    [0.00013075,-0.01651818,6.38550475],
    [4.95224596, -2.34187282,1.40722903],
    [-5.94120956,0.87524543,1.07437676],
    [-4.49880283,-3.12964705,1.41125120],
    [-0.75991343,-4.18114040,0.15674939],
    [2.20650584,-5.59236598,1.05533723],
    [0.96252338,4.36198028,-1.42717920],
    [4.81764990,3.19485539,-2.77622069],
    [-1.73138701,4.16471485,-3.46626689],
    [-2.55306625,2.91594904,4.51767596],
    [1.50206160,1.89218763,-4.09811156],
    [0.98329346,2.34564526,3.5974529],
    [4.46454220,-0.56369361,-3.47475802],
    [3.80849139,0.73094006,4.51682402],
    [-5.17525966,2.60395011,-2.77772471],
    [-2.37975873,0.36837952,-4.07698444],
    [-2.53045715,-0.33540920,3.59449188],
    [0.88144546,-2.23348348,-4.09931249],
    [1.54175224,-2.03600240,3.58132334],
    [3.31038688,-3.02228474,-1.45695989],
    [-4.27528734,-1.35796536,-1.44852313],
    [-2.72513737,-3.56991001,-3.48209886],
    [-1.25794100,-3.68787964,4.50492503],
    [0.35261163,-5.75579556,-2.79552608],
    [0.74123700,4.43697518,-4.02956973],
    [-0.64622546,4.44077927,3.74595547],
    [4.18920378,1.70873805,-4.66166552],
    [3.47773976,2.51701231,2.54165082],
    [-3.55947851,2.79106662,-4.65228042],
    [-3.92545305,1.73860005,2.54076867],
    [-0.00271263,0.01013464,-5.02518153],
    [-0.00332939,-0.00785361,2.48730355],
    [3.46776368,-2.83946171,-4.05736012],
    [4.17115409,-1.67896151,3.73089394],
    [-4.20196677,-1.57060486,-4.04933515],
    [-3.53060015,-2.79614465,3.73504438],
    [-0.62463325,-4.46372023,-4.67635054],
    [0.44748177,-4.28341329,2.52481152],
    [1.79334057,6.15314947,0.33208226],
    [5.55303831,3.65754425,-0.46937522],
    [-1.77455171,4.68869481,-0.74204629],
    [1.74568849,1.79694051,-1.26045912],
    [1.94632660,1.65173243,6.00004182],
    [4.95599592,-0.80678637,-0.75922862],
    [-5.94441654,2.98322195,-0.46759778],
    [-2.43853084,0.59877699,-1.25432294],
    [-2.41524712,0.83987083,5.99734552],
    [0.69721678,-2.42073954,-1.26854049],
    [0.46924459,-2.53096433,5.98775239],
    [4.43666586,-4.63906404,0.30520008],
    [-6.23050849,-1.53090073,0.31876024],
    [-3.17891818,-3.89146035,-0.76368896],
    [0.38866467,-6.64405085,-0.49626228]])
    print z.shape
 #   z = numpy.loadtxt('/home/chris/Desktop/q16_cd29/Cd29.xyz',skiprows = 2, usecols = (1,2,3), delimiter = ' ',unpack=False)
    chargemag = float(totalcharge)/charges
   
    print "number of charges",charges
    print 'total charge', totalcharge
    print 'magnitude of each charge', chargemag
    (x1,y1,z1) = tuple(numpy.mean(z,axis=0))
    center =  numpy.mean(z,axis=0)
    minradius = max(sqrt(sum((z[:]-center)**2,axis = 1)))
    radius = minradius + d
    chargecoordlist = pointsonsphere(charges,radius) + center

#    with open('/home/chris/Desktop/q16_cd29/Cd29.xyz','rb') as r:
#        xyzcoords = r.read()[3:]
#        r.close()
        
    totalcharge = str(totalcharge)
    print totalcharge
    newfolder = 'Q'+totalcharge+'d'+str(charges)+'r'+str(d)+'_Cd29S29'
    try:
        os.mkdir('/home/chris/Desktop/charges/'+newfolder)
    except:
        pass
    outputfilename = str('/home/chris/Desktop/charges/'+newfolder+'/'+newfolder+'.nw')
    newfolder = 'Q'+totalcharge+'d'+str(charges)+'r'+str(d)+'_Cd29S29'
    print "newfolder name", newfolder
    print "saving as", outputfilename
    with open(outputfilename,'wb') as f:
        f.write('''start "'''+newfolder+'''"
echo

memory total 4096 mb

charge 0
geometry
Cd 3.73479929 4.70148859 1.08453093
Cd -0.45777444 5.46389693 1.42922813
Cd 3.99256744 1.42776590 0.17101104
Cd -3.23253881 2.73557278 0.17389669
Cd 0.00277074 -0.00738694 -0.25202305
Cd 0.00013075 -0.01651818 6.38550475
Cd 4.95224596 -2.34187282 1.40722903
Cd -5.94120956 0.87524543 1.07437676
Cd -4.49880283 -3.12964705 1.41125120
Cd -0.75991343 -4.18114040 0.15674939
Cd 2.20650584 -5.59236598 1.05533723
Cd 0.96252338 4.36198028 -1.42717920
Cd 4.81764990 3.19485539 -2.77622069
Cd -1.73138701 4.16471485 -3.46626689
Cd -2.55306625 2.91594904 4.51767596
Cd 1.50206160 1.89218763 -4.09811156
Cd 0.98329346 2.34564526 3.59745292
Cd 4.46454220 -0.56369361 -3.47475802
Cd 3.80849139 0.73094006 4.51682402
Cd -5.17525966 2.60395011 -2.77772471
Cd -2.37975873 0.36837952 -4.07698444
Cd -2.53045715 -0.33540920 3.59449188
Cd 0.88144546 -2.23348348 -4.09931249
Cd 1.54175224 -2.03600240 3.58132334
Cd 3.31038688 -3.02228474 -1.45695989
Cd -4.27528734 -1.35796536 -1.44852313
Cd -2.72513737 -3.56991001 -3.48209886
Cd -1.25794100 -3.68787964 4.50492503
Cd 0.35261163 -5.75579556 -2.79552608
S 0.74123700 4.43697518 -4.02956973
S -0.64622546 4.44077927 3.74595547
S 4.18920378 1.70873805 -4.66166552
S 3.47773976 2.51701231 2.54165082
S -3.55947851 2.79106662 -4.65228042
S -3.92545305 1.73860005 2.54076867
S -0.00271263 0.01013464 -5.02518153
S -0.00332939 -0.00785361 2.48730355
S 3.46776368 -2.83946171 -4.05736012
S 4.17115409 -1.67896151 3.73089394
S -4.20196677 -1.57060486 -4.04933515
S -3.53060015 -2.79614465 3.73504438
S -0.62463325 -4.46372023 -4.67635054
S 0.44748177 -4.28341329 2.52481152
S 1.79334057 6.15314947 0.33208226
S 5.55303831 3.65754425 -0.46937522
S -1.77455171 4.68869481 -0.74204629
S 1.74568849 1.79694051 -1.26045912
S 1.94632660 1.65173243 6.00004182
S 4.95599592 -0.80678637 -0.75922862
S -5.94441654 2.98322195 -0.46759778
S -2.43853084 0.59877699 -1.25432294
S -2.41524712 0.83987083 5.99734552
S 0.69721678 -2.42073954 -1.26854049
S 0.46924459 -2.53096433 5.98775239
S 4.43666586 -4.63906404 0.30520008
S -6.23050849 -1.53090073 0.31876024
S -3.17891818 -3.89146035 -0.76368896
S 0.38866467 -6.64405085 -0.49626228
end 

bq\n''')
            
        for i in chargecoordlist:
                
               
                f.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+' ' +str(chargemag) + '\n')
                
        f.write('''
end

basis
Cd library "LANL2DZ ECP"
S library "LANL2DZ ECP"

end

ecp
   Cd library "LANL2DZ ECP"
   S library "LANL2DZ ECP"

end

dft
  xc b3lyp
  direct
  mult 1
  iterations 200
  convergence	energy 1E-5 	density 5E-6	gradient 5e-5 ncysh 200
end

dplot
  TITLE LUMO
  GAUSSIAN
  vectors ./'''+newfolder+'''.movecs
   LimitXYZ
 -10 10 50 
 -10 10 50
 -10  10  50
  orbitals view;1;262
  output lumo.cube
end
task tddft gradient
task dplot

dplot
  TITLE HOMO
  GAUSSIAN
  vectors ./'''+newfolder+'''.movecs
   LimitXYZ
 -10 10 50 
 -10 10 50
 -10  10  50
  orbitals view;1;261
  output homo.cube
end
task dplot''')
        f.close()
    if False: 
        fig=figure()
        ax = fig.add_subplot(111, projection='3d')
            
                
        Axes3D.scatter(ax,z[:,0],z[:,1],z[:,2],c='b')
        
        Axes3D.scatter(ax,chargecoordlist[:,0],chargecoordlist[:,1],chargecoordlist[:,2],c='r')
    return None

def CdSetddftcharged(charges,chargemag,totalcharge,d=0):  

    #### d distance of charges from surface in is in angstrom

    z = numpy.loadtxt('/home/chris/Desktop/CdSe29_finalgeometry.xyz',skiprows = 2, usecols = (1,2,3), delimiter = ' ',unpack=False)
    if z.shape[0]!=58:
        return 0
    chargemag = float(totalcharge)/charges
   
    print "number of charges",charges
    print 'total charge', totalcharge
    print 'magnitude of each charge', chargemag
    (x1,y1,z1) = tuple(numpy.mean(z,axis=0))
    center =  numpy.mean(z,axis=0)
    minradius = max(sqrt(sum((z[:]-center)**2,axis = 1)))
    radius = minradius + d
    chargecoordlist = pointsonsphere(charges,radius) + center

#    with open('/home/chris/Desktop/CdSe29_finalgeometry.xyz','rb') as r:
#        xyzcoords = r.read()[3:]
#        r.close()
        
    totalcharge = str(totalcharge)
    print totalcharge
    newfolder = 'tddftQ'+totalcharge+'d'+str(charges)+'r'+str(d)+'_Cd29Se29'
    try:
        os.mkdir('/home/chris/Desktop/charges/'+newfolder)
    except:
        pass
    outputfilename = str('/home/chris/Desktop/charges/'+newfolder+'/'+newfolder+'.nw')
    newfolder = 'tddftQ'+totalcharge+'d'+str(charges)+'r'+str(d)+'_Cd29Se29'
    print "newfolder name", newfolder
    print "saving as", outputfilename
    with open(outputfilename,'wb') as f:
        f.write('''start "'''+newfolder+'''"
echo

memory total 4096 mb

charge 0
geometry
symmetry c1\n''')
        for i in range(z.shape[0]):
            if i<29: f.write('Cd '+str(z[i,0])+' '+str(z[i,1])+' '+str(z[i,2])+'\n')
            elif i>30: f.write('Se '+str(z[i,0])+' '+str(z[i,1])+' '+str(z[i,2])+'\n')
        f.write('''end 

bq\n''')
            
        for i in chargecoordlist:
    
                f.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+' ' +str(chargemag) + '\n')
                
        f.write('''
end

basis
Cd library "LANL2DZ ECP"
Se library "LANL2DZ ECP"

end

ecp
   Cd library "LANL2DZ ECP"
   Se library "LANL2DZ ECP"

end

dft
  xc b3lyp
  direct
  mult 1
  iterations 200
  convergence	energy 1E-5 	density 5E-6	gradient 5e-5 ncysh 200
end

tddft
    nroots 10
    civecs
    notriplet

end

dplot
  TITLE LUMO
  GAUSSIAN
  vectors ./'''+newfolder+'''.movecs
   LimitXYZ
 -10 10 50 
 -10 10 50
 -10  10  50
  orbitals view;1;262
  output lumo.cube
end
task tddft energy
task dplot

dplot
  TITLE HOMO
  GAUSSIAN
  vectors ./'''+newfolder+'''movecs
   LimitXYZ
 -10 10 50 
 -10 10 50
 -10  10  50
  orbitals view;1;261
  output homo.cube
end
task dplot''')
        f.close()
    if False: 
        fig=figure()
        ax = fig.add_subplot(111, projection='3d')
            
                
        Axes3D.scatter(ax,z[:,0],z[:,1],z[:,2],c='b')
        
        Axes3D.scatter(ax,chargecoordlist[:,0],chargecoordlist[:,1],chargecoordlist[:,2],c='r')
    return None

def preptddftcharged(charges,chargemag,totalcharge,d=0):  

    #### d distance of charges from surface in is in angstrom

    z = numpy.loadtxt('/home/chris/Dropbox/DataWeiss/160113/OPT_LANL_B3LYP_nakedcluster/finalgeometry.xyz',skiprows = 2, usecols = (1,2,3), delimiter = ' ',unpack=False)
    chargemag = float(totalcharge)/charges
   
    print "number of charges",charges
    print 'total charge', totalcharge
    print 'magnitude of each charge', chargemag
    (x1,y1,z1) = tuple(numpy.mean(z,axis=0))
    center =  numpy.mean(z,axis=0)
    minradius = max(sqrt(sum((z[:]-center)**2,axis = 1)))
    radius = minradius + d
    chargecoordlist = pointsonsphere(charges,radius) + center

    
        
    totalcharge = str(totalcharge)
    print totalcharge
    newfolder = 'tddftQ'+totalcharge+'d'+str(charges)+'r'+str(d)+'_Cd29S29'
    try:
        os.mkdir('/home/chris/Desktop/charges/'+newfolder)
    except:
        pass
    outputfilename = str('/home/chris/Desktop/charges/'+newfolder+'/'+newfolder+'.nw')
    newfolder = 'tddftQ'+totalcharge+'d'+str(charges)+'r'+str(d)+'_Cd29S29'
    print "newfolder name", newfolder
    print "saving as", outputfilename
    with open(outputfilename,'wb') as f:
        f.write('''start "'''+newfolder+'''"
echo

memory total 4096 mb

charge 0
geometry
symmetry c1
Cd 3.73479929 4.70148859 1.08453093
Cd -0.45777444 5.46389693 1.42922813
Cd 3.99256744 1.42776590 0.17101104
Cd -3.23253881 2.73557278 0.17389669
Cd 0.00277074 -0.00738694 -0.25202305
Cd 0.00013075 -0.01651818 6.38550475
Cd 4.95224596 -2.34187282 1.40722903
Cd -5.94120956 0.87524543 1.07437676
Cd -4.49880283 -3.12964705 1.41125120
Cd -0.75991343 -4.18114040 0.15674939
Cd 2.20650584 -5.59236598 1.05533723
Cd 0.96252338 4.36198028 -1.42717920
Cd 4.81764990 3.19485539 -2.77622069
Cd -1.73138701 4.16471485 -3.46626689
Cd -2.55306625 2.91594904 4.51767596
Cd 1.50206160 1.89218763 -4.09811156
Cd 0.98329346 2.34564526 3.59745292
Cd 4.46454220 -0.56369361 -3.47475802
Cd 3.80849139 0.73094006 4.51682402
Cd -5.17525966 2.60395011 -2.77772471
Cd -2.37975873 0.36837952 -4.07698444
Cd -2.53045715 -0.33540920 3.59449188
Cd 0.88144546 -2.23348348 -4.09931249
Cd 1.54175224 -2.03600240 3.58132334
Cd 3.31038688 -3.02228474 -1.45695989
Cd -4.27528734 -1.35796536 -1.44852313
Cd -2.72513737 -3.56991001 -3.48209886
Cd -1.25794100 -3.68787964 4.50492503
Cd 0.35261163 -5.75579556 -2.79552608
S 0.74123700 4.43697518 -4.02956973
S -0.64622546 4.44077927 3.74595547
S 4.18920378 1.70873805 -4.66166552
S 3.47773976 2.51701231 2.54165082
S -3.55947851 2.79106662 -4.65228042
S -3.92545305 1.73860005 2.54076867
S -0.00271263 0.01013464 -5.02518153
S -0.00332939 -0.00785361 2.48730355
S 3.46776368 -2.83946171 -4.05736012
S 4.17115409 -1.67896151 3.73089394
S -4.20196677 -1.57060486 -4.04933515
S -3.53060015 -2.79614465 3.73504438
S -0.62463325 -4.46372023 -4.67635054
S 0.44748177 -4.28341329 2.52481152
S 1.79334057 6.15314947 0.33208226
S 5.55303831 3.65754425 -0.46937522
S -1.77455171 4.68869481 -0.74204629
S 1.74568849 1.79694051 -1.26045912
S 1.94632660 1.65173243 6.00004182
S 4.95599592 -0.80678637 -0.75922862
S -5.94441654 2.98322195 -0.46759778
S -2.43853084 0.59877699 -1.25432294
S -2.41524712 0.83987083 5.99734552
S 0.69721678 -2.42073954 -1.26854049
S 0.46924459 -2.53096433 5.98775239
S 4.43666586 -4.63906404 0.30520008
S -6.23050849 -1.53090073 0.31876024
S -3.17891818 -3.89146035 -0.76368896
S 0.38866467 -6.64405085 -0.49626228
end 

bq\n''')
            
        for i in chargecoordlist:
                
               
                f.write('Br '+str(i[0])+' '+str(i[1])+' '+str(i[2])+' \n')# +str(chargemag) + '\n')
                
        f.write('''
end

basis
Cd library "LANL2DZ ECP"
S library "LANL2DZ ECP"

end

ecp
   Cd library "LANL2DZ ECP"
   S library "LANL2DZ ECP"

end

dft
  xc b3lyp
  direct
  mult 1
  iterations 200
  convergence	energy 1E-5 	density 5E-6	gradient 5e-5 ncysh 200
end

tddft
    nroots 10
    civecs
    notriplet
    grad
        root 1
    end
end

dplot
  TITLE LUMO
  GAUSSIAN
  vectors ./'''+newfolder+'''.movecs
   LimitXYZ
 -10 10 50 
 -10 10 50
 -10  10  50
  orbitals view;1;262
  output lumo.cube
end
task tddft gradient
task dplot

dplot
  TITLE HOMO
  GAUSSIAN
  vectors ./'''+newfolder+'''movecs
   LimitXYZ
 -10 10 50 
 -10 10 50
 -10  10  50
  orbitals view;1;261
  output homo.cube
end
task dplot''')
        f.close()
    if False: 
        fig=figure()
        ax = fig.add_subplot(111, projection='3d')
            
                
        Axes3D.scatter(ax,z[:,0],z[:,1],z[:,2],c='b')
        
        Axes3D.scatter(ax,chargecoordlist[:,0],chargecoordlist[:,1],chargecoordlist[:,2],c='r')
    return None
    

def preparebatch():
#    prepdftcharged(1024,0,-1,d=4)
#    prepdftcharged(1024,0,-1,d=2)
#    prepdftcharged(1024,0,-1,d=1)
#    prepdftcharged(1024,0,-1,d=0.5)
#    prepdftcharged(256,0,-1,d=0.5)
#    prepdftcharged(16,0,-1,d=0.5)
#    prepdftcharged(8,0,-1,d=0.5)
#    prepdftcharged(4,0,-1,d=0.5)
#    prepdftcharged(1,0,-1,d=0.5)
#    prepdftcharged(1024,0,-1,d=0.5)
#    prepdftcharged(1024,0,-2,d=0.5)
#    prepdftcharged(1024,0,-3,d=0.5)
#    prepdftcharged(1024,0,-4,d=0.5)
#    prepdftcharged(1024,0,-5,d=0.5)
#    prepdftcharged(1024,0,-6,d=0.5)
#    prepdftcharged(1024,0,-8,d=0.5)
#    prepdftcharged(1024,0,-12,d=0.5)
#    prepdftcharged(1024,0,-16,d=0.5)
#    prepdftcharged(1024,0,-1,d=1)
#    prepdftcharged(1024,0,-2,d=1)
#    prepdftcharged(1024,0,-3,d=1)
#    prepdftcharged(1024,0,-4,d=1)
#    prepdftcharged(1024,0,-5,d=1)
#    prepdftcharged(1024,0,-6,d=1)
#    prepdftcharged(1024,0,-8,d=1)
#    prepdftcharged(1024,0,-12,d=1)
#    prepdftcharged(1024,0,-16,d=1)
#    prepdftcharged(512,0,-1,d=1)
#    prepdftcharged(256,0,-1,d=1)
#    prepdftcharged(128,0,-1,d=1)
#    prepdftcharged(64,0,-1,d=1)
#    prepdftcharged(32,0,-1,d=1)
#    prepdftcharged(16,0,-1,d=1)
#    prepdftcharged(8,0,-1,d=1)
#    prepdftcharged(4,0,-1,d=1)
#    prepdftcharged(2,0,-1,d=1)
#    prepdftcharged(1,0,-1,d=1)
#    preptddftcharged(1024,0,-1,d=1)
#    preptddftcharged(1024,0,-4,d=1)
#    preptddftcharged(1024,0,-2,d=1)
#    preptddftcharged(1024,0,-16,d=1)
#    preptddftcharged(1024,0,-8,d=1)  
#
#    prepdftcharged(1024,0,-16,d=2)
#    prepdftcharged(1024,0,-8,d=2)
#    prepdftcharged(1024,0,-4,d=2)
#    prepdftcharged(1024,0,-2,d=2)
#    prepdftcharged(1024,0,-1,d=2)
#    prepdftcharged(1024,0,-256,d=4)
#    prepdftcharged(1024,0,-128,d=4)
#    prepdftcharged(1024,0,-64,d=4)
#    prepdftcharged(1024,0,-32,d=4)
#    prepdftcharged(1024,0,-16,d=4)

    

#    prepdftcharged(1024,0,1,d=1)
#    prepdftcharged(1024,0,2,d=1)
#    prepdftcharged(1024,0,3,d=1)
#    prepdftcharged(1024,0,4,d=1)
#    prepdftcharged(1024,0,5,d=1)
#    prepdftcharged(1024,0,6,d=1)
#    prepdftcharged(1024,0,8,d=1)
#    prepdftcharged(1024,0,12,d=1)
#    prepdftcharged(1024,0,16,d=1)

    CdSetddftcharged(1024,0,0,d=1)
    CdSetddftcharged(1024,0,-1,d=1)
    CdSetddftcharged(1024,0,-2,d=1)
    CdSetddftcharged(1024,0,-4,d=1)
    CdSetddftcharged(1024,0,-8,d=1)
    CdSetddftcharged(1024,0,1,d=1)
    CdSetddftcharged(1024,0,2,d=1)
    CdSetddftcharged(1024,0,4,d=1)
    CdSetddftcharged(1024,0,8,d=1)

    
    return 0
    


#x = array([[11.26 ,21.1026656797 ,17.5],
#[13.34,17.5,17.5],
#[15.42,21.1026656797,17.5],
#[15.42,13.8973343203,17.5],
#[17.5,17.5,17.5],
#[17.5,17.5,24.256],
#[19.58,21.1026656797,17.5],
#[17.5,10.2946686405,17.5],
#[19.58,13.8973343203,17.5],
#[21.66,17.5,17.5],
#[23.74,21.1026656797,17.5],
#[11.24544,18.6996876714,14.122],
#[13.32544,22.3023533511,14.122],
#[13.32544,15.0970219916,14.122],
#[13.32544,15.0970219916,20.878],
#[15.40544,18.6996876714,14.122],
#[15.40544,18.6996876714,20.878],
#[17.48544,22.3023533511,14.122],
#[17.48544,22.3023533511,20.878],
#[15.40544,11.4943563119,14.122],
#[17.48544,15.0970219916,14.122],
#[17.48544,15.0970219916,20.878],
#[19.56544,18.6996876714,14.122],
#[19.56544,18.6996876714,20.878],
#[21.64544,22.3023533511,14.122],
#[19.56544,11.4943563119,14.122],
#[21.64544,15.0970219916,14.122],
#[21.64544,15.0970219916,20.878],
#[23.72544,18.6996876714,14.122],
#[13.34,17.5,13.2775],
#[13.34,17.5,20.0335],
#[15.42,21.1026656797,13.2775],
#[15.42,21.1026656797,20.033],
#[15.42,13.8973343203,13.2775],
#[15.42,13.8973343203,20.0335],
#[17.5,17.5,13.2775],
#[17.5,17.5,20.0335],
#[19.58,21.1026656797,13.2775],
#[19.58,21.1026656797,20.0335],
#[19.58,13.8973343203,13.2775],
#[19.58,13.8973343203,20.0335],
#[21.66,17.5,13.2775],
#[21.66,17.5,20.0335],
#[11.259792,18.7007684711,16.6555],
#[13.339792,22.3034341508,16.6555],
#[13.339792,15.0981027913,16.6555],
#[15.419792,18.7007684711,16.6555],
#[15.419792,18.7007684711,23.4115],
#[17.499792,22.3034341508,16.6555],
#[15.419792,11.4954371116,16.6555],
#[17.499792,15.0981027913,16.6555],
#[17.499792,15.0981027913,23.4115],
#[19.579792,18.7007684711,16.6555],
#[19.579792,18.7007684711,23.4115],
#[21.659792,22.3034341508,16.6555],
#[19.579792,11.4954371116,16.6555],
#[21.659792,15.0981027913,16.6555],
#[23.739792,18.7007684711,16.6555]])
#
#x*=[[1,1,7.0152/6.756]]
#x*=[[4.2985/4.16,1,1]]
#f = open('/home/chris/Desktop/CdSeParticle2.xyz','wb')  
#f.write('58\nCdSe\n')  
#for i in x[0:29]:
#    f.write('Cd '+str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
#for i in x[29:]:
#    f.write('Se '+str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
#    
#f.close()

