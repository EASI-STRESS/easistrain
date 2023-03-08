import os
import sys
import time
import configparser
from datetime import datetime
from .func_integration_2D import integration_2D

# time0 = time.time()

# if (len(sys.argv) != 11 or sys.argv[1] == 'help'):
# 	print('#### Arguments to be given: 1)root // 2)h5file // 3)scan // 4)detector_name // 5)path of poni file // 6)npt_rad // 7)npt_azim // 8)unit(2th_deg or 2th_rad or q_A^-1...) //9)dark image if exist if not give 0//10)mask image if exist if not give 0')
# else:
# 	integration_2D(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])

# time1 = time.time()

# print('total time: ' + str(time1-time0) + ' seconds')


# -----------------------------delete all above except import when the script below is finished


# Use: python exe_integration_2D.py filename.ini

time0 = time.time()  # Record the time when starting the script

config = configparser.ConfigParser()  #
config.read(
    sys.argv[1]
)  # read of the first argument given in the terminal (the configuration file 'filename.ini')

root_data = config.get("arguments", "rootDir")
h5file = config.get("arguments", "h5File")
scan = config.get("arguments", "scan")

if (
    config.get("arguments", "numScanStart") == "None"
    or config.get("arguments", "numScanEnd") == "None"
):
    numScan = None
else:
    numScan = (
        int(config.get("arguments", "numScanStart")),
        int(config.get("arguments", "numScanEnd")),
    )

detector_name = config.get("arguments", "detectorName")
poni_file = config.get("arguments", "poniFile")
npt_rad = config.get("arguments", "nptRad")
npt_azim = config.get("arguments", "nptAzim")
x_unit = config.get("arguments", "xUnit")

if config.get("arguments", "imDark") == "None":
    im_dark = None
else:
    im_dark = config.get("arguments", "imDark")

if config.get("arguments", "imMask") == "None":
    im_mask = None
else:
    im_mask = config.get("arguments", "imMask")

if (
    config.get("arguments", "radRangeStart") == "None"
    or config.get("arguments", "radRangeEnd") == "None"
):
    rad_range = None
else:
    rad_range = (
        float(config.get("arguments", "radRangeStart")),
        float(config.get("arguments", "radRangeEnd")),
    )

if (
    config.get("arguments", "azimRangeStart") == "None"
    or config.get("arguments", "azimRangeEnd") == "None"
):
    azim_range = None
else:
    azim_range = (
        float(config.get("arguments", "azimRangeStart")),
        float(config.get("arguments", "azimRangeEnd")),
    )

if config.get("arguments", "errorModel") == "None":
    errorModel = None
else:
    errorModel = config.get("arguments", "errorModel")

if config.get("arguments", "imFlat") == "None":
    imFlat = None
else:
    imFlat = config.get("arguments", "imFlat")

if config.get("arguments", "gon1") == "None":
    gon1 = None
else:
    gon1 = config.get("arguments", "gon1")

if config.get("arguments", "gon2") == "None":
    gon2 = None
else:
    gon2 = config.get("arguments", "gon2")

if config.get("arguments", "gon3") == "None":
    gon3 = None
else:
    gon3 = config.get("arguments", "gon3")

if config.get("arguments", "chiGon1") == "None":
    chiGon1 = None
else:
    chiGon1 = config.get("arguments", "chiGon1")

if config.get("arguments", "omegaGon2") == "None":
    omegaGon2 = None
else:
    omegaGon2 = config.get("arguments", "omegaGon2")

if config.get("arguments", "phiGon3") == "None":
    phiGon3 = None
else:
    phiGon3 = config.get("arguments", "phiGon3")

print(azim_range, rad_range)
print(type(azim_range), type(rad_range))

if os.path.dirname(root_data):
    os.makedirs(os.path.dirname(root_data), exist_ok=True)

with open(os.path.join(root_data, "exe_integration.log"), "w") as fwlog:
    fwlog.write("INTEGRATION 2D LOG FILE\n")
    fwlog.write("Date and time : " + str(datetime.now()) + "\n")

with open(os.path.join(root_data, "exe_integration.log"), "a") as falog:
    falog.write("#$#$#$#$#$#$#The arguments used for the integration are below: \n")
    falog.write("root_data = " + root_data + "\n")
    falog.write("h5file = " + h5file + "\n")
    falog.write("scan = " + scan + "\n")
    falog.write("(numScanStart,numScanEnd) = (" + str(numScan) + ")\n")
    falog.write("detector_name = " + detector_name + "\n")
    falog.write("poni_file = " + poni_file + "\n")
    falog.write("npt_rad = " + npt_rad + "\n")
    falog.write("npt_azim = " + npt_azim + "\n")
    falog.write("x_unit = " + x_unit + "\n")
    falog.write("im_dark = " + str(im_dark) + "\n")
    falog.write("im_mask = " + str(im_mask) + "\n")
    falog.write("(rad_rangeStart,rad_rangeEnd) = (" + str(rad_range) + "\n")
    falog.write("(azim_rangeStart,azim_rangeEnd) = (" + str(azim_range) + "\n")
    falog.write("errorModel = " + str(errorModel) + "\n")
    falog.write("imFlat = " + str(imFlat) + "\n")
    falog.write("gon1 = " + str(gon1) + "\n")
    falog.write("gon2 = " + str(gon2) + "\n")
    falog.write("gon3 = " + str(gon3) + "\n")
    falog.write("chiGon1 = " + str(chiGon1) + "\n")
    falog.write("omegaGon2 = " + str(omegaGon2) + "\n")
    falog.write("phiGon3 = " + str(phiGon3) + "\n")
    falog.write("************____________________________________**************\n")

    integration_2D(
        root_data,
        h5file,
        scan,
        numScan,
        detector_name,
        poni_file,
        npt_rad,
        npt_azim,
        x_unit,
        im_dark,
        im_mask,
        rad_range,
        azim_range,
        errorModel,
        imFlat,
        gon1,
        gon2,
        gon3,
        chiGon1,
        omegaGon2,
        phiGon3,
    )
    time1 = time.time()
    print("total time: " + str(time1 - time0) + " seconds")
    falog.write("total time: " + str(time1 - time0) + " seconds\n")
