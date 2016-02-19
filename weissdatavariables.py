# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 20:00:42 2016

@author: chris
"""

import os, sys


__datafolder = '/home/chris/Dropbox/DataWeiss'
__scriptsfolder = '/home/chris/Dropbox/PyScripts/weiss'







if __scriptsfolder not in sys.path:
    
    sys.path.append(__scriptsfolder)

os.chdir(__datafolder)
         
    