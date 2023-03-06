# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 18:45:26 2021

@author: slim
"""

import numpy

# ----definition of the normalized x, y and z vectors ----
vx = numpy.array([1, 0, 0])
vy = numpy.array([0, 1, 0])
vz = numpy.array([0, 0, 1])


# ----Rotation matrix around x, x is in the direction of the beam ----
def matrotx(rx):
    rotx = numpy.array(
        [
            [1, 0, 0],
            [0, numpy.cos(rx), -numpy.sin(rx)],
            [0, numpy.sin(rx), numpy.cos(rx)],
        ]
    )
    return rotx


# ----Rotation matrix around y, y is perpendicular to the beam and in its plane ----
def matroty(ry):
    roty = numpy.array(
        [
            [numpy.cos(ry), 0, numpy.sin(ry)],
            [0, 1, 0],
            [-numpy.sin(ry), 0, numpy.cos(ry)],
        ]
    )
    return roty


# ----Rotation matrix around z, z iz perpendicular to the beam and out of the beam plane ----
def matrotz(rz):
    rotz = numpy.array(
        [
            [numpy.cos(rz), -numpy.sin(rz), 0],
            [numpy.sin(rz), numpy.cos(rz), 0],
            [0, 0, 1],
        ]
    )
    return rotz


# ----Rotation matrix for 3 rotations: first around x, second around y and third around z  ----
def matrotxyz(rx, ry, rz):
    rotxyz = numpy.dot(matrotz(rz), numpy.dot(matroty(ry), matrotx(rx)))
    return rotxyz


# ----Norm of a vector ----
def normvec(cx, cy, cz):
    normv = numpy.sqrt(pow(cx, 2) + pow(cy, 2) + pow(cz, 2))
    return normv


# ----angle de diffraction theta ----
def tth(cx, cy, cz, rx, ry, rz):
    theta = 0.5 * numpy.arctan2(
        numpy.sqrt(
            pow(
                numpy.dot(
                    numpy.dot(matrotxyz(rx, ry, rz), numpy.array([cx, cy, cz])), vy
                ),
                2,
            )
            + pow(
                numpy.dot(
                    numpy.dot(matrotxyz(rx, ry, rz), numpy.array([cx, cy, cz])), vz
                ),
                2,
            )
        ),
        numpy.dot(numpy.dot(matrotxyz(rx, ry, rz), numpy.array([cx, cy, cz])), vx),
    )
    return theta


def azim(cx, cy, cz, rx, ry, rz):
    azimuth = numpy.arctan2(
        numpy.dot(numpy.dot(matrotxyz(rx, ry, rz), numpy.array([cx, cy, cz])), vz),
        numpy.dot(numpy.dot(matrotxyz(rx, ry, rz), numpy.array([cx, cy, cz])), vy),
    )
    return azimuth
