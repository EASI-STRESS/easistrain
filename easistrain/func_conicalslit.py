# -*- coding: utf-8 -*-
"""
Created on Thu May 20 17:32:56 2021

@author: slim
"""
import numpy as np


# This function determine the distance between the centre of the silt and the conical aperture
# lss is the distance between the slit and the sample
# r is the radius (distance between the centre of the silt and the aperture)

def slitradius(a,lss,tth):
    r = lss*np.tan(2*tth)
    return r
    
# This function determine the depth resolution of the conical slit (length of the gauge volume)
# dpr is the depth resolution (length of the gauge volume)
# lsd is the distance between the slit and the detector
# lss is the distance between the slit and the sample
# dss is the distance travelled by the diffraction cone before imtersectiong the slit
# dsd is the distance travelled by the diffraction between the slit and the detector
# dslit is the opening of the slit
# dfoc is the defocus of the beam
# dpsd is the spatial resolution of the detector (pixel size)
# The formula was extracted from the paper of Lienet et al [Lienert,Mat.Res.Soc.Symp.Proc,vol.590,2000]
def lengthgv(tth,lss,lsd,dslit,dfoc,dpsd):
    dss = lss/np.cos(2*tth)
    dsd = ((lss+lsd)/np.cos(2*tth))-dss
    dpr = (1/np.tan(2*tth))*np.sqrt(((((dss+dsd)*dslit)/dsd)**2)+(dfoc**2)+(((dss*dpsd)/dsd)**2))
    return dpr

# This function determine the depth resolution of the conical slit (length of the gauge volume)
# lss is the distance between the slit and the sample
# rs is the radius of the slit (distance between the centre of the slit and the opening)
# bs is the beam size in the veritical direction
# dslit is the opening of the slit
# I write this formula (It is an approximation if the Lienert formula)
def gvlength(rs,dslit,bs,lss):
    gvl=((dslit+bs)*lss)/(rs+(dslit/2))
    return gvl
    