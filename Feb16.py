# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 15:26:51 2015

@author: chris

"""
import pandas
        
def EDSproc(a = None):
    
    silver = array([])  
    os.chdir('/home/chris/Documents/DataWeiss/150216_EDS/Chris thompson')
    if a == None:
        a=load('/home/chris/Documents/DataWeiss/150216_EDS/Chris thompson/4-datacube.npy')
    else:  
        a = copy(a)

    a = a.reshape((2048,256,-1))
    
    silver = sum(a[610:620,:,:],axis=0)
    figure()
    title('silver')
    plt.imshow(silver,cmap = cm.gist_heat_r)
    savefig('7-silver.png')

    selenium = sum(a[300:314,:,:],axis=0)
    figure()
    title('selenium')
    plt.imshow(selenium, cmap=cm.gist_heat_r)
    savefig('7-selenium.png')
    
    
#    subplot(223)
#    hist(reshape(silver,(-1,)),bins =[-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5])# [-0.5,0.5,100,200,300,400,500,600,700,800,900,1000])
#    ylim(0,10)
    
    figure()
    title('spectrum')
    plot(linspace(-0.1,10.135,2048),sum(a, axis=(1,2)))
    
        
    return 0
    

def imtoarray(filename):
    os.chdir('/home/chris/Documents/DataWeiss/150216_EDS/Chris thompson')
    from PIL import Image
    im = Image.open(filename)
    
    a=array(im.getdata()).reshape(im.size+(3,))  ### 3 for RGB image?
    return a
    
def trythis():    
    myfile ='/home/chris/Documents/DataWeiss/150216_EDS/Chris thompson/4-smartmap.raw'
    
    with open(myfile) as f:
        f.seek(0)
        a = fromiter((int(i) for i in binascii.b2a_hex(f.read())),dtype=int)
        
        f.close()
        a=respace(a(2048,256,-1))
        save('/home/chris/Documents/DataWeiss/150216_EDS/Chris thompson/4-smartmap.npy',a)
    return 0
    
def re(x):
    z = '/home/chris/Documents/DataWeiss/150316'
    os.chdir(z)    
    for x in os.listdir(z):
        if '.txt' in x and 'notes' not in x:
            f = open(z+'/'+x,'rb')
            b = f.read()
            b = b.replace('\t',',')
            f.close()
            f=open(z+'/'+x,'wb')
            f.write(b)
            f.close()
            
            r = loadtxt(z+'/'+x,dtype = float,delimiter = ',')# [',','\t',' ']
            
            r[:,0] = 10**7/(r[:,0]- 10**7/785)  
           
            savetxt(z+'/'+x,r,delimiter = ',')
    print os.listdir('.')
    return 0
    
def re2(x):
    z = '/home/chris/Documents/DataWeiss/150318'
    if '.txt' in x and 'notes' not in x:
       
        r = loadtxt(z+'/'+x,dtype = float,delimiter = ',')# [',','\t',' ']
          
        r[:,0] = 10**7/r[:,0] - 10**7/785
        savetxt(z+'/'+x,r,delimiter = ',')
    plot(r[:,0],r[:,1])
    print os.listdir(z)
    return 0
    