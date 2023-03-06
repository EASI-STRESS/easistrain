#!/usr/bin/python3
"""
usage : python imgshow.py image_file.edf [max_int_value]
"""

import pylab
import fabio
import numpy
import sys

img = fabio.open(sys.argv[1])
pixels = img.data
maxpixels = numpy.amax(pixels)
minpixels = numpy.amin(pixels)
print(maxpixels, minpixels)
if numpy.size(sys.argv) == 4:
    pylab.imshow(
        img.data,
        vmin=int(sys.argv[2]),
        vmax=int(sys.argv[3]),
        cmap="gray_r",
        origin="lower",
    )
if numpy.size(sys.argv) == 3:
    pylab.imshow(img.data, vmax=int(sys.argv[2]), cmap="gray_r", origin="lower")
if numpy.size(sys.argv) == 2:
    pylab.imshow(img.data, vmax=10000, cmap="gray_r", origin="lower")
cb = pylab.colorbar()
pylab.show()
