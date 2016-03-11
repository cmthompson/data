# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 20:00:42 2016

@author: chris

Weissdatavariables.py contains important stuff for using the data workup scripts I wrote.  It should be imported by any 
of the data processing scripts.
"""

import os, sys


__datafolder = '/home/chris/Dropbox/DataWeiss'
__scriptsfolder = '/home/chris/Dropbox/PyScripts/weiss'







if __scriptsfolder not in sys.path:
    
    sys.path.append(__scriptsfolder)

os.chdir(__datafolder)
         
    