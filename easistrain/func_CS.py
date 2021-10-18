# -*- coding: utf-8 -*-
"""
Created on Wed May 26 17:29:11 2021

@author: slim
"""

import numpy as np
from tthdspacing import *
from conicalslit import *
from plot import *

# Units of the different quantities
### energy = keV
### distance = micron
### angle =rad
###

energy = np.linspace(50, 150, 10000)  # energy (keV) of x-rays

# The function 'bccslit' is for a BCC structure (just put the desired lattice parameter, ap)

# Definition of the different quantities needed for calculation
# ap = 0.0002855 # lattice parameter of the desired BCC structure
# sddistance = 1000000 # slit to detector distance, sddistance
# slitopening = 25 # the slit opening
# beamsize = 50 # the size of the used beam
# pixelsize = 200 # the pixel size of the used detector


def bccslit(ap, sddistance, slitopening, beamsize, pixelsize, title):
    # Calculation of d-spacing and bragg angle theta
    d110, tth110 = cubicdspacing(energy, ap, 1, 1, 0)
    d200, tth200 = cubicdspacing(energy, ap, 2, 0, 0)
    d211, tth211 = cubicdspacing(energy, ap, 2, 1, 1)
    d220, tth220 = cubicdspacing(energy, ap, 2, 2, 0)
    d310, tth310 = cubicdspacing(energy, ap, 3, 1, 0)
    d222, tth222 = cubicdspacing(energy, ap, 2, 2, 2)
    d321, tth321 = cubicdspacing(energy, ap, 3, 2, 1)
    d400, tth400 = cubicdspacing(energy, ap, 4, 0, 0)
    ### radius of the slit at the 110 peak ###
    sltr11050 = slitradius(ap, 50000, tth110)
    sltr11060 = slitradius(ap, 60000, tth110)
    sltr11070 = slitradius(ap, 70000, tth110)
    sltr11080 = slitradius(ap, 80000, tth110)
    sltr11090 = slitradius(ap, 90000, tth110)
    sltr110100 = slitradius(ap, 100000, tth110)
    ### radius of the slit at the 200 peak ###
    sltr20050 = slitradius(ap, 50000, tth200)
    sltr20060 = slitradius(ap, 60000, tth200)
    sltr20070 = slitradius(ap, 70000, tth200)
    sltr20080 = slitradius(ap, 80000, tth200)
    sltr20090 = slitradius(ap, 90000, tth200)
    sltr200100 = slitradius(ap, 100000, tth200)
    ### radius of the slit at the 211 peak ###
    sltr21150 = slitradius(ap, 50000, tth211)
    sltr21160 = slitradius(ap, 60000, tth211)
    sltr21170 = slitradius(ap, 70000, tth211)
    sltr21180 = slitradius(ap, 80000, tth211)
    sltr21190 = slitradius(ap, 90000, tth211)
    sltr211100 = slitradius(ap, 100000, tth211)
    ### radius of the slit at the 220 peak ###
    sltr22050 = slitradius(ap, 50000, tth220)
    sltr22060 = slitradius(ap, 60000, tth220)
    sltr22070 = slitradius(ap, 70000, tth220)
    sltr22080 = slitradius(ap, 80000, tth220)
    sltr22090 = slitradius(ap, 90000, tth220)
    sltr220100 = slitradius(ap, 100000, tth220)
    ### radius of the slit at the 310 peak ###
    sltr31050 = slitradius(ap, 50000, tth310)
    sltr31060 = slitradius(ap, 60000, tth310)
    sltr31070 = slitradius(ap, 70000, tth310)
    sltr31080 = slitradius(ap, 80000, tth310)
    sltr31090 = slitradius(ap, 90000, tth310)
    sltr310100 = slitradius(ap, 100000, tth310)
    ### radius of the slit at the 222 peak ###
    sltr22250 = slitradius(ap, 50000, tth222)
    sltr22260 = slitradius(ap, 60000, tth222)
    sltr22270 = slitradius(ap, 70000, tth222)
    sltr22280 = slitradius(ap, 80000, tth222)
    sltr22290 = slitradius(ap, 90000, tth222)
    sltr222100 = slitradius(ap, 100000, tth222)
    ### radius of the slit at the 321 peak ###
    sltr32150 = slitradius(ap, 50000, tth321)
    sltr32160 = slitradius(ap, 60000, tth321)
    sltr32170 = slitradius(ap, 70000, tth321)
    sltr32180 = slitradius(ap, 80000, tth321)
    sltr32190 = slitradius(ap, 90000, tth321)
    sltr321100 = slitradius(ap, 100000, tth321)
    ### radius of the slit at the 400 peak ###
    sltr40050 = slitradius(ap, 50000, tth400)
    sltr40060 = slitradius(ap, 60000, tth400)
    sltr40070 = slitradius(ap, 70000, tth400)
    sltr40080 = slitradius(ap, 80000, tth400)
    sltr40090 = slitradius(ap, 90000, tth400)
    sltr400100 = slitradius(ap, 100000, tth400)
    ### Plotting of the radius for the 110 peak ###
    horaxis01 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis01 = np.array(
        [
            0.001 * sltr11050,
            0.001 * sltr11060,
            0.001 * sltr11070,
            0.001 * sltr11080,
            0.001 * sltr11090,
            0.001 * sltr110100,
        ]
    )
    legends01 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis01,
        veraxis01,
        "energy (keV)",
        "slitradius110 (mm)",
        "-",
        legends01,
        "radius_BCC (110)",
        title + " (110)",
    )
    ### Plotting of the radius for the 200 peak ###
    horaxis02 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis02 = np.array(
        [
            0.001 * sltr20050,
            0.001 * sltr20060,
            0.001 * sltr20070,
            0.001 * sltr20080,
            0.001 * sltr20090,
            0.001 * sltr200100,
        ]
    )
    legends02 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis02,
        veraxis02,
        "energy (keV)",
        "slitradius200 (mm)",
        "-",
        legends02,
        "radius_BCC (200)",
        title + " (200)",
    )
    ### Plotting of the radius for the 211 peak ###
    horaxis03 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis03 = np.array(
        [
            0.001 * sltr21150,
            0.001 * sltr21160,
            0.001 * sltr21170,
            0.001 * sltr21180,
            0.001 * sltr21190,
            0.001 * sltr211100,
        ]
    )
    legends03 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis03,
        veraxis03,
        "energy (keV)",
        "slitradius211 (mm)",
        "-",
        legends03,
        "radius_BCC (211)",
        title + " (211)",
    )
    ### Plotting of the radius for the 220 peak ###
    horaxis04 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis04 = np.array(
        [
            0.001 * sltr22050,
            0.001 * sltr22060,
            0.001 * sltr22070,
            0.001 * sltr22080,
            0.001 * sltr22090,
            0.001 * sltr220100,
        ]
    )
    legends04 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis04,
        veraxis04,
        "energy (keV)",
        "slitradius220 (mm)",
        "-",
        legends04,
        "radius_BCC (220)",
        title + " (220)",
    )
    ### Plotting of the radius for the 310 peak ###
    horaxis05 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis05 = np.array(
        [
            0.001 * sltr31050,
            0.001 * sltr31060,
            0.001 * sltr31070,
            0.001 * sltr31080,
            0.001 * sltr31090,
            0.001 * sltr310100,
        ]
    )
    legends05 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis05,
        veraxis05,
        "energy (keV)",
        "slitradius310 (mm)",
        "-",
        legends05,
        "radius_BCC (310)",
        title + " (310)",
    )
    ### Plotting of the radius for the 222 peak ###
    horaxis06 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis06 = np.array(
        [
            0.001 * sltr22250,
            0.001 * sltr22260,
            0.001 * sltr22270,
            0.001 * sltr22280,
            0.001 * sltr22290,
            0.001 * sltr222100,
        ]
    )
    legends06 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis06,
        veraxis06,
        "energy (keV)",
        "slitradius222 (mm)",
        "-",
        legends06,
        "radius_BCC (222)",
        title + " (222)",
    )
    ### Plotting of the radius for the 321 peak ###
    horaxis07 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis07 = np.array(
        [
            0.001 * sltr32150,
            0.001 * sltr32160,
            0.001 * sltr32170,
            0.001 * sltr32180,
            0.001 * sltr32190,
            0.001 * sltr321100,
        ]
    )
    legends07 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis07,
        veraxis07,
        "energy (keV)",
        "slitradius321 (mm)",
        "-",
        legends07,
        "radius_BCC (321)",
        title + " (321)",
    )
    ### Plotting of the radius for the 400 peak ###
    horaxis08 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis08 = np.array(
        [
            0.001 * sltr40050,
            0.001 * sltr40060,
            0.001 * sltr40070,
            0.001 * sltr40080,
            0.001 * sltr40090,
            0.001 * sltr400100,
        ]
    )
    legends08 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis08,
        veraxis08,
        "energy (keV)",
        "slitradius211 (mm)",
        "-",
        legends08,
        "radius_BCC (400)",
        title + " (400)",
    )
    ## Gauge volume length of the 110 peak ####
    gvl11050 = lengthgv(tth110, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl11060 = lengthgv(tth110, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl11070 = lengthgv(tth110, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl11080 = lengthgv(tth110, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl11090 = lengthgv(tth110, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl110100 = lengthgv(tth110, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Gauge volume length of the 200 peak ####
    gvl20050 = lengthgv(tth200, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl20060 = lengthgv(tth200, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl20070 = lengthgv(tth200, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl20080 = lengthgv(tth200, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl20090 = lengthgv(tth200, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl200100 = lengthgv(tth200, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Gauge volume length of the 211 peak ####
    gvl21150 = lengthgv(tth211, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl21160 = lengthgv(tth211, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl21170 = lengthgv(tth211, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl21180 = lengthgv(tth211, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl21190 = lengthgv(tth211, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl211100 = lengthgv(tth211, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Gauge volume length of the 220 peak ####
    gvl22050 = lengthgv(tth220, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl22060 = lengthgv(tth220, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl22070 = lengthgv(tth220, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl22080 = lengthgv(tth220, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl22090 = lengthgv(tth220, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl220100 = lengthgv(tth220, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Gauge volume length of the 310 peak ####
    gvl31050 = lengthgv(tth310, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl31060 = lengthgv(tth310, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl31070 = lengthgv(tth310, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl31080 = lengthgv(tth310, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl31090 = lengthgv(tth310, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl310100 = lengthgv(tth310, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Gauge volume length of the 222 peak ####
    gvl22250 = lengthgv(tth222, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl22260 = lengthgv(tth222, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl22270 = lengthgv(tth222, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl22280 = lengthgv(tth222, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl22290 = lengthgv(tth222, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl222100 = lengthgv(tth222, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Gauge volume length of the 321 peak ####
    gvl32150 = lengthgv(tth321, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl32160 = lengthgv(tth321, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl32170 = lengthgv(tth321, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl32180 = lengthgv(tth321, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl32190 = lengthgv(tth321, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl321100 = lengthgv(tth321, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Gauge volume length of the 400 peak ####
    gvl40050 = lengthgv(tth400, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl40060 = lengthgv(tth400, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl40070 = lengthgv(tth400, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl40080 = lengthgv(tth400, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl40090 = lengthgv(tth400, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl400100 = lengthgv(tth400, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Plotting of the gauge volume length for the 110 peak ###
    horaxis11 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis11 = np.array([gvl11050, gvl11060, gvl11070, gvl11080, gvl11090, gvl110100])
    legends11 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis11,
        veraxis11,
        "energy (keV)",
        "gauge volume length 110 (μm)",
        "-",
        legends11,
        "gvl_BCC (110)",
        title + " (110)",
    )
    ### Plotting of the gauge volume length for the 200 peak ###
    horaxis12 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis12 = np.array([gvl20050, gvl20060, gvl20070, gvl20080, gvl20090, gvl200100])
    legends12 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis12,
        veraxis12,
        "energy (keV)",
        "gauge volume length 200 (μm)",
        "-",
        legends12,
        "gvl_BCC (200)",
        title + " (200)",
    )
    ### Plotting of the gauge volume length for the 211 peak ###
    horaxis13 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis13 = np.array([gvl21150, gvl21160, gvl21170, gvl21180, gvl21190, gvl211100])
    legends13 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis13,
        veraxis13,
        "energy (keV)",
        "gauge volume length 211 (μm)",
        "-",
        legends13,
        "gvl_BCC (211)",
        title + "(211)",
    )
    ### Plotting of the gauge volume length for the 220 peak ###
    horaxis14 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis14 = np.array([gvl22050, gvl22060, gvl22070, gvl22080, gvl22090, gvl220100])
    legends14 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis14,
        veraxis14,
        "energy (keV)",
        "gauge volume length 220 (μm)",
        "-",
        legends14,
        "gvl_BCC (220)",
        title + " (220)",
    )
    ### Plotting of the gauge volume length for the 310 peak ###
    horaxis15 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis15 = np.array([gvl31050, gvl31060, gvl31070, gvl31080, gvl31090, gvl310100])
    legends15 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis15,
        veraxis15,
        "energy (keV)",
        "gauge volume length 310 (μm)",
        "-",
        legends15,
        "gvl_BCC (310)",
        title + " (310)",
    )
    ### Plotting of the gauge volume length for the 222 peak ###
    horaxis16 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis16 = np.array([gvl22250, gvl22260, gvl22270, gvl22280, gvl22290, gvl222100])
    legends16 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis16,
        veraxis16,
        "energy (keV)",
        "gauge volume length 222 (μm)",
        "-",
        legends16,
        "gvl_BCC (222)",
        title + "(222)",
    )
    ### Plotting of the gauge volume length for the 321 peak ###
    horaxis17 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis17 = np.array([gvl32150, gvl32160, gvl32170, gvl32180, gvl32190, gvl321100])
    legends17 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis17,
        veraxis17,
        "energy (keV)",
        "gauge volume length 321 (μm)",
        "-",
        legends17,
        "gvl_BCC (321)",
        title + " (321)",
    )
    ### Plotting of the gauge volume length for the 400 peak ###
    horaxis18 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis18 = np.array([gvl40050, gvl40060, gvl40070, gvl40080, gvl40090, gvl400100])
    legends18 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis18,
        veraxis18,
        "energy (keV)",
        "gauge volume length 400 (μm)",
        "-",
        legends18,
        "gvl_BCC (400)",
        title + " (400)",
    )
    ### Plotting of the radius of all planes for one slit to sample distance
    horaxis09 = np.array(
        [energy, energy, energy, energy, energy, energy, energy, energy]
    )
    veraxis09 = np.array(
        [
            0.001 * sltr110100,
            0.001 * sltr200100,
            0.001 * sltr211100,
            0.001 * sltr220100,
            0.001 * sltr310100,
            0.001 * sltr222100,
            0.001 * sltr321100,
            0.001 * sltr400100,
        ]
    )
    legends09 = np.array(
        ["(110)", "(200)", "(211)", "(220)", "(310)", "(222)", "(321)", "(400)"]
    )
    showplot(
        horaxis09,
        veraxis09,
        "energy (keV)",
        "slit radius (mm)",
        "-",
        legends09,
        "radius_BCC_lss100mm",
        title,
    )
    ### Plotting of the length of the gauge volume of all planes for one slit to sample distance
    horaxis19 = np.array(
        [energy, energy, energy, energy, energy, energy, energy, energy]
    )
    veraxis19 = np.array(
        [
            gvl110100,
            gvl200100,
            gvl211100,
            gvl220100,
            gvl310100,
            gvl222100,
            gvl321100,
            gvl400100,
        ]
    )
    legends19 = np.array(
        ["(110)", "(200)", "(211)", "(220)", "(310)", "(222)", "(321)", "(400)"]
    )
    showplot(
        horaxis19,
        veraxis19,
        "energy (keV)",
        "gauge volume length (μm)",
        "-",
        legends19,
        "gvl_BCC_lss100mm",
        title,
    )
    return


### fcc structure
def fccslit(ap, sddistance, slitopening, beamsize, pixelsize, title):
    # Calculation of d-spacing and bragg angle theta
    d111, tth111 = cubicdspacing(energy, ap, 1, 1, 1)
    d200, tth200 = cubicdspacing(energy, ap, 2, 0, 0)
    d220, tth220 = cubicdspacing(energy, ap, 2, 2, 0)
    d311, tth311 = cubicdspacing(energy, ap, 3, 1, 1)
    d222, tth222 = cubicdspacing(energy, ap, 2, 2, 2)
    d400, tth400 = cubicdspacing(energy, ap, 4, 0, 0)
    d331, tth331 = cubicdspacing(energy, ap, 3, 3, 1)
    d420, tth420 = cubicdspacing(energy, ap, 4, 2, 0)
    d422, tth422 = cubicdspacing(energy, ap, 4, 2, 2)
    ### radius of the slit at the 111 peak ###
    sltr11150 = slitradius(ap, 50000, tth111)
    sltr11160 = slitradius(ap, 60000, tth111)
    sltr11170 = slitradius(ap, 70000, tth111)
    sltr11180 = slitradius(ap, 80000, tth111)
    sltr11190 = slitradius(ap, 90000, tth111)
    sltr111100 = slitradius(ap, 100000, tth111)
    ### radius of the slit at the 200 peak ###
    sltr20050 = slitradius(ap, 50000, tth200)
    sltr20060 = slitradius(ap, 60000, tth200)
    sltr20070 = slitradius(ap, 70000, tth200)
    sltr20080 = slitradius(ap, 80000, tth200)
    sltr20090 = slitradius(ap, 90000, tth200)
    sltr200100 = slitradius(ap, 100000, tth200)
    ### radius of the slit at the 220 peak ###
    sltr22050 = slitradius(ap, 50000, tth220)
    sltr22060 = slitradius(ap, 60000, tth220)
    sltr22070 = slitradius(ap, 70000, tth220)
    sltr22080 = slitradius(ap, 80000, tth220)
    sltr22090 = slitradius(ap, 90000, tth220)
    sltr220100 = slitradius(ap, 100000, tth220)
    ### radius of the slit at the 311 peak ###
    sltr31150 = slitradius(ap, 50000, tth311)
    sltr31160 = slitradius(ap, 60000, tth311)
    sltr31170 = slitradius(ap, 70000, tth311)
    sltr31180 = slitradius(ap, 80000, tth311)
    sltr31190 = slitradius(ap, 90000, tth311)
    sltr311100 = slitradius(ap, 100000, tth311)
    ### radius of the slit at the 222 peak ###
    sltr22250 = slitradius(ap, 50000, tth222)
    sltr22260 = slitradius(ap, 60000, tth222)
    sltr22270 = slitradius(ap, 70000, tth222)
    sltr22280 = slitradius(ap, 80000, tth222)
    sltr22290 = slitradius(ap, 90000, tth222)
    sltr222100 = slitradius(ap, 100000, tth222)
    ### radius of the slit at the 311 peak ###
    sltr40050 = slitradius(ap, 50000, tth400)
    sltr40060 = slitradius(ap, 60000, tth400)
    sltr40070 = slitradius(ap, 70000, tth400)
    sltr40080 = slitradius(ap, 80000, tth400)
    sltr40090 = slitradius(ap, 90000, tth400)
    sltr400100 = slitradius(ap, 100000, tth400)
    ### radius of the slit at the 331 peak ###
    sltr33150 = slitradius(ap, 50000, tth331)
    sltr33160 = slitradius(ap, 60000, tth331)
    sltr33170 = slitradius(ap, 70000, tth331)
    sltr33180 = slitradius(ap, 80000, tth331)
    sltr33190 = slitradius(ap, 90000, tth331)
    sltr331100 = slitradius(ap, 100000, tth331)
    ### radius of the slit at the 420 peak ###
    sltr42050 = slitradius(ap, 50000, tth420)
    sltr42060 = slitradius(ap, 60000, tth420)
    sltr42070 = slitradius(ap, 70000, tth420)
    sltr42080 = slitradius(ap, 80000, tth420)
    sltr42090 = slitradius(ap, 90000, tth420)
    sltr420100 = slitradius(ap, 100000, tth420)
    ### radius of the slit at the 422 peak ###
    sltr42250 = slitradius(ap, 50000, tth422)
    sltr42260 = slitradius(ap, 60000, tth422)
    sltr42270 = slitradius(ap, 70000, tth422)
    sltr42280 = slitradius(ap, 80000, tth422)
    sltr42290 = slitradius(ap, 90000, tth422)
    sltr422100 = slitradius(ap, 100000, tth422)
    ### Plotting of the radius for the 111 peak ###
    horaxis21 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis21 = np.array(
        [
            0.001 * sltr11150,
            0.001 * sltr11160,
            0.001 * sltr11170,
            0.001 * sltr11180,
            0.001 * sltr11190,
            0.001 * sltr111100,
        ]
    )
    legends21 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis21,
        veraxis21,
        "energy (keV)",
        "slitradius111 (mm)",
        "-",
        legends21,
        "radius_FCC (111)",
        title + " (111)",
    )
    ### Plotting of the radius for the 200 peak ###
    horaxis22 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis22 = np.array(
        [
            0.001 * sltr20050,
            0.001 * sltr20060,
            0.001 * sltr20070,
            0.001 * sltr20080,
            0.001 * sltr20090,
            0.001 * sltr200100,
        ]
    )
    legends22 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis22,
        veraxis22,
        "energy (keV)",
        "slitradius200 (mm)",
        "-",
        legends22,
        "radius_FCC (200)",
        title + " (200)",
    )
    ### Plotting of the radius for the 220 peak ###
    horaxis23 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis23 = np.array(
        [
            0.001 * sltr22050,
            0.001 * sltr22060,
            0.001 * sltr22070,
            0.001 * sltr22080,
            0.001 * sltr22090,
            0.001 * sltr220100,
        ]
    )
    legends23 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis23,
        veraxis23,
        "energy (keV)",
        "slitradius220 (mm)",
        "-",
        legends23,
        "radius_FCC (220)",
        title + " (220)",
    )
    ### Plotting of the radius for the 311 peak ###
    horaxis24 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis24 = np.array(
        [
            0.001 * sltr31150,
            0.001 * sltr31160,
            0.001 * sltr31170,
            0.001 * sltr31180,
            0.001 * sltr31190,
            0.001 * sltr311100,
        ]
    )
    legends24 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis24,
        veraxis24,
        "energy (keV)",
        "slitradius311 (mm)",
        "-",
        legends24,
        "radius_FCC (311)",
        title + "(311)",
    )
    ### Plotting of the radius for the 222 peak ###
    horaxis25 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis25 = np.array(
        [
            0.001 * sltr22250,
            0.001 * sltr22260,
            0.001 * sltr22270,
            0.001 * sltr22280,
            0.001 * sltr22290,
            0.001 * sltr222100,
        ]
    )
    legends25 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis25,
        veraxis25,
        "energy (keV)",
        "slitradius222 (mm)",
        "-",
        legends25,
        "radius_FCC (222)",
        title + " (222)",
    )
    ### Plotting of the radius for the 400 peak ###
    horaxis26 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis26 = np.array(
        [
            0.001 * sltr40050,
            0.001 * sltr40060,
            0.001 * sltr40070,
            0.001 * sltr40080,
            0.001 * sltr40090,
            0.001 * sltr400100,
        ]
    )
    legends26 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis26,
        veraxis26,
        "energy (keV)",
        "slitradius400 (mm)",
        "-",
        legends26,
        "radius_FCC (400)",
        title + " (400)",
    )
    ### Plotting of the radius for the 331 peak ###
    horaxis27 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis27 = np.array(
        [
            0.001 * sltr33150,
            0.001 * sltr33160,
            0.001 * sltr33170,
            0.001 * sltr33180,
            0.001 * sltr33190,
            0.001 * sltr331100,
        ]
    )
    legends27 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis27,
        veraxis27,
        "energy (keV)",
        "slitradius331 (mm)",
        "-",
        legends27,
        "radius_FCC (331)",
        title + " (331)",
    )
    ### Plotting of the radius for the 420 peak ###
    horaxis28 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis28 = np.array(
        [
            0.001 * sltr42050,
            0.001 * sltr42060,
            0.001 * sltr42070,
            0.001 * sltr42080,
            0.001 * sltr42090,
            0.001 * sltr420100,
        ]
    )
    legends28 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis28,
        veraxis28,
        "energy (keV)",
        "slitradius420 (mm)",
        "-",
        legends28,
        "radius_FCC (420)",
        title + " (420)",
    )
    ### Plotting of the radius for the 422 peak ###
    horaxis29 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis29 = np.array(
        [
            0.001 * sltr42250,
            0.001 * sltr42260,
            0.001 * sltr42270,
            0.001 * sltr42280,
            0.001 * sltr42290,
            0.001 * sltr422100,
        ]
    )
    legends29 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis29,
        veraxis29,
        "energy (keV)",
        "slitradius422 (mm)",
        "-",
        legends29,
        "radius_FCC (422)",
        title + " (422)",
    )
    ## Gauge volume length of the 111 peak ####
    gvl11150 = lengthgv(tth111, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl11160 = lengthgv(tth111, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl11170 = lengthgv(tth111, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl11180 = lengthgv(tth111, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl11190 = lengthgv(tth111, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl111100 = lengthgv(tth111, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 200 peak ####
    gvl20050 = lengthgv(tth200, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl20060 = lengthgv(tth200, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl20070 = lengthgv(tth200, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl20080 = lengthgv(tth200, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl20090 = lengthgv(tth200, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl200100 = lengthgv(tth200, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 220 peak ####
    gvl22050 = lengthgv(tth220, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl22060 = lengthgv(tth220, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl22070 = lengthgv(tth220, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl22080 = lengthgv(tth220, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl22090 = lengthgv(tth220, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl220100 = lengthgv(tth220, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 311 peak ####
    gvl31150 = lengthgv(tth311, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl31160 = lengthgv(tth311, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl31170 = lengthgv(tth311, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl31180 = lengthgv(tth311, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl31190 = lengthgv(tth311, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl311100 = lengthgv(tth311, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 222 peak ####
    gvl22250 = lengthgv(tth222, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl22260 = lengthgv(tth222, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl22270 = lengthgv(tth222, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl22280 = lengthgv(tth222, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl22290 = lengthgv(tth222, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl222100 = lengthgv(tth222, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 400 peak ####
    gvl40050 = lengthgv(tth400, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl40060 = lengthgv(tth400, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl40070 = lengthgv(tth400, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl40080 = lengthgv(tth400, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl40090 = lengthgv(tth400, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl400100 = lengthgv(tth400, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 331 peak ####
    gvl33150 = lengthgv(tth331, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl33160 = lengthgv(tth331, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl33170 = lengthgv(tth331, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl33180 = lengthgv(tth331, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl33190 = lengthgv(tth331, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl331100 = lengthgv(tth331, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 420 peak ####
    gvl42050 = lengthgv(tth420, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl42060 = lengthgv(tth420, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl42070 = lengthgv(tth420, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl42080 = lengthgv(tth420, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl42090 = lengthgv(tth420, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl420100 = lengthgv(tth420, 100000, sddistance, slitopening, beamsize, pixelsize)
    ## Gauge volume length of the 422 peak ####
    gvl42250 = lengthgv(tth422, 50000, sddistance, slitopening, beamsize, pixelsize)
    gvl42260 = lengthgv(tth422, 60000, sddistance, slitopening, beamsize, pixelsize)
    gvl42270 = lengthgv(tth422, 70000, sddistance, slitopening, beamsize, pixelsize)
    gvl42280 = lengthgv(tth422, 80000, sddistance, slitopening, beamsize, pixelsize)
    gvl42290 = lengthgv(tth422, 90000, sddistance, slitopening, beamsize, pixelsize)
    gvl422100 = lengthgv(tth422, 100000, sddistance, slitopening, beamsize, pixelsize)
    ### Plotting of the gauge volume length for the 111 peak ###
    horaxis31 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis31 = np.array([gvl11150, gvl11160, gvl11170, gvl11180, gvl11190, gvl111100])
    legends31 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis31,
        veraxis31,
        "energy (keV)",
        "gauge volume length 111 (μm)",
        "-",
        legends31,
        "gvl_FCC (111)",
        title + " (111)",
    )
    ### Plotting of the gauge volume length for the 200 peak ###
    horaxis32 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis32 = np.array([gvl20050, gvl20060, gvl20070, gvl20080, gvl20090, gvl200100])
    legends32 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis32,
        veraxis32,
        "energy (keV)",
        "gauge volume length 200 (μm)",
        "-",
        legends32,
        "gvl_FCC (200)",
        title + " (200)",
    )
    ### Plotting of the gauge volume length for the 220 peak ###
    horaxis33 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis33 = np.array([gvl22050, gvl22060, gvl22070, gvl22080, gvl22090, gvl220100])
    legends33 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis33,
        veraxis33,
        "energy (keV)",
        "gauge volume length 220 (μm)",
        "-",
        legends33,
        "gvl_FCC (220)",
        title + " (220)",
    )
    ### Plotting of the gauge volume length for the 311 peak ###
    horaxis34 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis34 = np.array([gvl31150, gvl31160, gvl31170, gvl31180, gvl31190, gvl311100])
    legends34 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis34,
        veraxis34,
        "energy (keV)",
        "gauge volume length 311 (μm)",
        "-",
        legends34,
        "gvl_FCC (311)",
        title + " (311)",
    )
    ### Plotting of the gauge volume length for the 222 peak ###
    horaxis35 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis35 = np.array([gvl22250, gvl22260, gvl22270, gvl22280, gvl22290, gvl222100])
    legends35 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis35,
        veraxis35,
        "energy (keV)",
        "gauge volume length 222 (μm)",
        "-",
        legends35,
        "gvl_FCC (222)",
        title + " (222)",
    )
    ### Plotting of the gauge volume length for the 400 peak ###
    horaxis36 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis36 = np.array([gvl40050, gvl40060, gvl40070, gvl40080, gvl40090, gvl400100])
    legends36 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis36,
        veraxis36,
        "energy (keV)",
        "gauge volume length 400 (μm)",
        "-",
        legends36,
        "gvl_FCC (400)",
        title + " (400)",
    )
    ### Plotting of the gauge volume length for the 331 peak ###
    horaxis37 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis37 = np.array([gvl33150, gvl33160, gvl33170, gvl33180, gvl33190, gvl331100])
    legends37 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis37,
        veraxis37,
        "energy (keV)",
        "gauge volume length 331 (μm)",
        "-",
        legends37,
        "gvl_FCC (331)",
        title + " (331)",
    )
    ### Plotting of the gauge volume length for the 420 peak ###
    horaxis38 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis38 = np.array([gvl42050, gvl42060, gvl42070, gvl42080, gvl42090, gvl420100])
    legends38 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis38,
        veraxis38,
        "energy (keV)",
        "gauge volume length 420 (μm)",
        "-",
        legends38,
        "gvl_FCC (420)",
        title + " (420)",
    )
    ### Plotting of the gauge volume length for the 422 peak ###
    horaxis39 = np.array([energy, energy, energy, energy, energy, energy])
    veraxis39 = np.array([gvl42250, gvl42260, gvl42270, gvl42280, gvl42290, gvl422100])
    legends39 = np.array(
        [
            "lss = 50 mm",
            "lss = 60 mm",
            "lss = 70 mm",
            "lss = 80 mm",
            "lss = 90 mm",
            "lss = 100 mm",
        ]
    )
    showplot(
        horaxis39,
        veraxis39,
        "energy (keV)",
        "gauge volume length 422 (μm)",
        "-",
        legends39,
        "gvl_FCC (422)",
        title + " (422)",
    )
    ### Plotting of the radius of all planes for one slit to sample distance
    horaxis029 = np.array(
        [energy, energy, energy, energy, energy, energy, energy, energy, energy]
    )
    veraxis029 = np.array(
        [
            0.001 * sltr111100,
            0.001 * sltr200100,
            0.001 * sltr220100,
            0.001 * sltr311100,
            0.001 * sltr222100,
            0.001 * sltr400100,
            0.001 * sltr331100,
            0.001 * sltr420100,
            0.001 * sltr422100,
        ]
    )
    legends029 = np.array(
        [
            "(111)",
            "(200)",
            "(220)",
            "(311)",
            "(222)",
            "(400)",
            "(331)",
            "(420)",
            "(422)",
        ]
    )
    showplot(
        horaxis029,
        veraxis029,
        "energy (keV)",
        "slitradius(mm)",
        "-",
        legends029,
        "radius_FCC_lss100mm",
        title,
    )
    ### Plotting of the gauge volume length of all planes for one slit to sample distance
    horaxis029 = np.array(
        [energy, energy, energy, energy, energy, energy, energy, energy, energy]
    )
    veraxis029 = np.array(
        [
            gvl111100,
            gvl200100,
            gvl220100,
            gvl311100,
            gvl222100,
            gvl400100,
            gvl331100,
            gvl420100,
            gvl422100,
        ]
    )
    legends029 = np.array(
        [
            "(111)",
            "(200)",
            "(220)",
            "(311)",
            "(222)",
            "(400)",
            "(331)",
            "(420)",
            "(422)",
        ]
    )
    showplot(
        horaxis029,
        veraxis029,
        "energy (keV)",
        "gauge volume length (μm)",
        "-",
        legends029,
        "gvl_FCC_lss100mm",
        title,
    )
    return


#### This function calculate the radius of the slit for a hkl plane
def slitradiushkl(energy, ap, h, k, l, lss):
    # Calculation of d-spacing and bragg angle theta
    d, tth = cubicdspacing(energy, ap, h, k, l)
    ### radius of the slit at the 111 peak ###
    sltrhkl = slitradius(ap, lss, tth)
    return sltrhkl


#### This function plot the radius of two structures BCC and FCC
def csFCCBCC(aBCC, aFCC, phasebcc, phasefcc, lss, title):
    ### BCC structure ####
    dbcc110, tthbcc110 = cubicdspacing(energy, aBCC, 1, 1, 0)
    dbcc200, tthbcc200 = cubicdspacing(energy, aBCC, 2, 0, 0)
    dbcc211, tthbcc211 = cubicdspacing(energy, aBCC, 2, 1, 1)
    dbcc220, tthbcc220 = cubicdspacing(energy, aBCC, 2, 2, 0)
    dbcc310, tthbcc310 = cubicdspacing(energy, aBCC, 3, 1, 0)
    dbcc222, tthbcc222 = cubicdspacing(energy, aBCC, 2, 2, 2)
    dbcc321, tthbcc321 = cubicdspacing(energy, aBCC, 3, 2, 1)
    dbcc400, tthbcc400 = cubicdspacing(energy, aBCC, 4, 0, 0)
    bcc110 = slitradius(aBCC, lss, tthbcc110)
    bcc200 = slitradius(aBCC, lss, tthbcc200)
    bcc211 = slitradius(aBCC, lss, tthbcc211)
    bcc220 = slitradius(aBCC, lss, tthbcc220)
    bcc310 = slitradius(aBCC, lss, tthbcc310)
    bcc222 = slitradius(aBCC, lss, tthbcc222)
    bcc321 = slitradius(aBCC, lss, tthbcc321)
    bcc400 = slitradius(aBCC, lss, tthbcc400)
    plt.figure(num=title, figsize=(10, 8))
    plt.plot(energy, 0.001 * bcc110, "-", label=phasebcc + "(110)", linewidth=3)
    plt.plot(energy, 0.001 * bcc200, "-", label=phasebcc + "(200)", linewidth=3)
    plt.plot(energy, 0.001 * bcc211, "-", label=phasebcc + "(211)", linewidth=3)
    plt.plot(energy, 0.001 * bcc220, "-", label=phasebcc + "(220)", linewidth=3)
    plt.plot(energy, 0.001 * bcc310, "-", label=phasebcc + "(310)", linewidth=3)
    plt.plot(energy, 0.001 * bcc222, "-", label=phasebcc + "(222)", linewidth=3)
    plt.plot(energy, 0.001 * bcc321, "-", label=phasebcc + "(321)", linewidth=3)
    plt.plot(energy, 0.001 * bcc400, "-", label=phasebcc + "(400)", linewidth=3)
    ### FCC structure ###
    dfcc111, tthfcc111 = cubicdspacing(energy, aFCC, 1, 1, 1)
    dfcc200, tthfcc200 = cubicdspacing(energy, aFCC, 2, 0, 0)
    dfcc220, tthfcc220 = cubicdspacing(energy, aFCC, 2, 2, 0)
    dfcc311, tthfcc311 = cubicdspacing(energy, aFCC, 3, 1, 1)
    dfcc222, tthfcc222 = cubicdspacing(energy, aFCC, 2, 2, 2)
    dfcc400, tthfcc400 = cubicdspacing(energy, aFCC, 4, 0, 0)
    dfcc331, tthfcc331 = cubicdspacing(energy, aFCC, 3, 3, 1)
    dfcc420, tthfcc420 = cubicdspacing(energy, aFCC, 4, 2, 0)
    dfcc422, tthfcc422 = cubicdspacing(energy, aFCC, 4, 2, 2)
    fcc111 = slitradius(aFCC, lss, tthfcc111)
    fcc200 = slitradius(aFCC, lss, tthfcc200)
    fcc220 = slitradius(aFCC, lss, tthfcc220)
    fcc311 = slitradius(aFCC, lss, tthfcc311)
    fcc222 = slitradius(aFCC, lss, tthfcc222)
    fcc400 = slitradius(aFCC, lss, tthfcc400)
    fcc331 = slitradius(aFCC, lss, tthfcc331)
    fcc420 = slitradius(aFCC, lss, tthfcc420)
    fcc422 = slitradius(aFCC, lss, tthfcc422)
    plt.plot(energy, 0.001 * fcc111, "-.", label=phasefcc + "(111)", linewidth=3)
    plt.plot(energy, 0.001 * fcc200, "-.", label=phasefcc + "(200)", linewidth=3)
    plt.plot(energy, 0.001 * fcc220, "-.", label=phasefcc + "(220)", linewidth=3)
    plt.plot(energy, 0.001 * fcc311, "-.", label=phasefcc + "(311)", linewidth=3)
    plt.plot(energy, 0.001 * fcc222, "-.", label=phasefcc + "(222)", linewidth=3)
    plt.plot(energy, 0.001 * fcc400, "-.", label=phasefcc + "(400)", linewidth=3)
    plt.plot(energy, 0.001 * fcc331, "-.", label=phasefcc + "(331)", linewidth=3)
    plt.plot(energy, 0.001 * fcc420, "-.", label=phasefcc + "(420)", linewidth=3)
    plt.plot(energy, 0.001 * fcc422, "-.", label=phasefcc + "(422)", linewidth=3)
    ### PLOT ###
    plt.xlabel("energy (keV)", family="sans-serif", fontsize=28)
    plt.ylabel("slit radius (mm)", family="sans-serif", fontsize=28)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(loc="best", fontsize=22)
    plt.grid()
    plt.title(title, fontsize=30)
    plt.show()
    # plt.savefig(title,dpi=200)
    plt.close()
    return


#### This function plot the radius of two FCC structures
def csFCCFCC(a1, a2, phase1, phase2, lss, title):
    ### BCC structure ####
    d1fcc111, tth1fcc111 = cubicdspacing(energy, a1, 1, 1, 1)
    d1fcc200, tth1fcc200 = cubicdspacing(energy, a1, 2, 0, 0)
    d1fcc220, tth1fcc220 = cubicdspacing(energy, a1, 2, 2, 0)
    d1fcc311, tth1fcc311 = cubicdspacing(energy, a1, 3, 1, 1)
    d1fcc222, tth1fcc222 = cubicdspacing(energy, a1, 2, 2, 2)
    d1fcc400, tth1fcc400 = cubicdspacing(energy, a1, 4, 0, 0)
    d1fcc331, tth1fcc331 = cubicdspacing(energy, a1, 3, 3, 1)
    d1fcc420, tth1fcc420 = cubicdspacing(energy, a1, 4, 2, 0)
    d1fcc422, tth1fcc422 = cubicdspacing(energy, a1, 4, 2, 2)
    onefcc111 = slitradius(a1, lss, tth1fcc111)
    onefcc200 = slitradius(a1, lss, tth1fcc200)
    onefcc220 = slitradius(a1, lss, tth1fcc220)
    onefcc311 = slitradius(a1, lss, tth1fcc311)
    onefcc222 = slitradius(a1, lss, tth1fcc222)
    onefcc400 = slitradius(a1, lss, tth1fcc400)
    onefcc331 = slitradius(a1, lss, tth1fcc331)
    onefcc420 = slitradius(a1, lss, tth1fcc420)
    onefcc422 = slitradius(a1, lss, tth1fcc422)
    plt.figure(num=title, figsize=(10, 8))
    plt.plot(energy, 0.001 * onefcc111, "-", label=phase1 + "(111)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc200, "-", label=phase1 + "(200)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc220, "-", label=phase1 + "(220)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc311, "-", label=phase1 + "(311)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc222, "-", label=phase1 + "(222)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc400, "-", label=phase1 + "(400)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc331, "-", label=phase1 + "(331)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc420, "-", label=phase1 + "(420)", linewidth=3)
    plt.plot(energy, 0.001 * onefcc422, "-", label=phase1 + "(422)", linewidth=3)
    ### FCC structure ###
    d2fcc111, tth2fcc111 = cubicdspacing(energy, a2, 1, 1, 1)
    d2fcc200, tth2fcc200 = cubicdspacing(energy, a2, 2, 0, 0)
    d2fcc220, tth2fcc220 = cubicdspacing(energy, a2, 2, 2, 0)
    d2fcc311, tth2fcc311 = cubicdspacing(energy, a2, 3, 1, 1)
    d2fcc222, tth2fcc222 = cubicdspacing(energy, a2, 2, 2, 2)
    d2fcc400, tth2fcc400 = cubicdspacing(energy, a2, 4, 0, 0)
    d2fcc331, tth2fcc331 = cubicdspacing(energy, a2, 3, 3, 1)
    d2fcc420, tth2fcc420 = cubicdspacing(energy, a2, 4, 2, 0)
    d2fcc422, tth2fcc422 = cubicdspacing(energy, a2, 4, 2, 2)
    twofcc111 = slitradius(a2, lss, tth2fcc111)
    twofcc200 = slitradius(a2, lss, tth2fcc200)
    twofcc220 = slitradius(a2, lss, tth2fcc220)
    twofcc311 = slitradius(a2, lss, tth2fcc311)
    twofcc222 = slitradius(a2, lss, tth2fcc222)
    twofcc400 = slitradius(a2, lss, tth2fcc400)
    twofcc331 = slitradius(a2, lss, tth2fcc331)
    twofcc420 = slitradius(a2, lss, tth2fcc420)
    twofcc422 = slitradius(a2, lss, tth2fcc422)
    plt.plot(energy, 0.001 * twofcc111, "-.", label=phase2 + "(111)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc200, "-.", label=phase2 + "(200)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc220, "-.", label=phase2 + "(220)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc311, "-.", label=phase2 + "(311)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc222, "-.", label=phase2 + "(222)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc400, "-.", label=phase2 + "(400)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc331, "-.", label=phase2 + "(331)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc420, "-.", label=phase2 + "(420)", linewidth=3)
    plt.plot(energy, 0.001 * twofcc422, "-.", label=phase2 + "(422)", linewidth=3)
    ### PLOT ###
    plt.xlabel("energy (keV)", family="sans-serif", fontsize=28)
    plt.ylabel("slit radius (mm)", family="sans-serif", fontsize=28)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(loc="best", fontsize=22)
    plt.grid()
    plt.title(title, fontsize=30)
    plt.show()
    # plt.savefig(title,dpi=200)
    plt.close()
    return
