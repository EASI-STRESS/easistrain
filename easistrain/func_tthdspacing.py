# -*- coding: utf-8 -*-
"""
Created on Thu May 20 09:24:46 2021

@author: slim
"""

import numpy as np

# This function calculte the d-spacing of a plane of a cubic structure
# a is the lattice parameter of the structure in micron
# h, k and l are the miller indices of the crystallographic plane
# e is the eneregy in keV

def cubicdspacing(e,a,h,k,l):
    d = a/(np.sqrt((h**2)+(k**2)+(l**2)))
    theta = np.arcsin(0.0012398/(2*e*d))
    return d, theta

# This function calculate the d-spacing of a plane of a hexagonal structure
# a and c are the lattice parameters of the structure in micron
# h, k and l are the miller indices of the crystallographic plane
# e is the eneregy in keV

def hexdspacing(e,a,c,h,k,l):
    d = np.sqrt(1/((4/3)*(((h**2)+h*k+(k**2))/(a**2))+((l**2)/(c**2))))
    theta = np.arcsin(0.0012398/(2*e*d))
    return d, theta