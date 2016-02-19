# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 21:09:51 2016

@author: chris
"""




import sys, inspect

#def is_mod_function(mod, func):
#    return inspect.isfunction(func) and inspect.getmodule(func) == mod
#
#def list_functions(mod):
#    return [func.__name__ for func in mod.__dict__.itervalues() 
#            if is_mod_function(mod, func)]
#
#
#print 'functions in current module:\n', list_functions(sys.modules[__name__])
#print 'functions in inspect module:\n', list_functions(inspect)

import inspect
import imp

os.chdir('/home/chris/Dropbox/PyScripts/weiss/')
with open('/home/chris/Desktop/filelist.txt','wb') as s:
    for afile in os.listdir('.'):
        
        if afile[-5:] =='ES.py':
        
            continue
        elif afile[-3:]=='.py':
    
            
            
            foo = imp.load_source('some','/home/chris/Dropbox/PyScripts/weiss/'+afile)
            s.write('\n'+afile[:-3])
            
            for name, data in inspect.getmembers(foo):
                if name == '__builtins__':
                    continue
                if inspect.getmodule(data)==foo:
                    
                    s.write('\n\t'+str(data.__name__)+':'+str(data.__doc__))
            print 'imported', afile     
           
            
    s.close()           
        
            
                
           
            