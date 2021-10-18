# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 18:45:26 2021

@author: slim
"""

import numpy as np

### definition of the normalized x, y and z vectors ###
vx = np.array([1, 0, 0])
vy = np.array([0, 1, 0])
vz = np.array([0, 0, 1])
### Rotation matrix around x, x is in the direction of the beam ###
def matrotx(rx):
    rotx = np.array(
        [[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]]
    )
    return rotx


### Rotation matrix around y, y is perpendicular to the beam and in its plane ###
def matroty(ry):
    roty = np.array(
        [[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]]
    )
    return roty


### Rotation matrix around z, z iz perpendicular to the beam and out of the beam plane ###
def matrotz(rz):
    rotz = np.array(
        [[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]]
    )
    return rotz


### Rotation matrix for 3 rotations: first around x, second around y and third around z  ###
def matrotxyz(rx, ry, rz):
    rotxyz = np.dot(matrotz(rz), np.dot(matroty(ry), matrotx(rx)))
    return rotxyz


### Norm of a vector ###
def normvec(cx, cy, cz):
    normv = np.sqrt(pow(cx, 2) + pow(cy, 2) + pow(cz, 2))
    return normv


### angle de diffraction theta ###
def tth(cx, cy, cz, rx, ry, rz):
    theta = 0.5 * np.arctan2(
        np.sqrt(
            pow(np.dot(np.dot(matrotxyz(rx, ry, rz), np.array([cx, cy, cz])), vy), 2)
            + pow(np.dot(np.dot(matrotxyz(rx, ry, rz), np.array([cx, cy, cz])), vz), 2)
        ),
        np.dot(np.dot(matrotxyz(rx, ry, rz), np.array([cx, cy, cz])), vx),
    )
    return theta


def azim(cx, cy, cz, rx, ry, rz):
    azimuth = np.arctan2(
        np.dot(np.dot(matrotxyz(rx, ry, rz), np.array([cx, cy, cz])), vz),
        np.dot(np.dot(matrotxyz(rx, ry, rz), np.array([cx, cy, cz])), vy),
    )
    return azimuth
