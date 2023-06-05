Unordered notes
===============

Task 1: `calibrationEDD`
++++++++++++++++++++++++

Purpose: calibration of the conversion from channels to energy with X-ray source

.. math::

    \text{E} = p_0 x^2 + p_1 x + p_2

Procedure:

- fit peaks in channels
- fit energies (e.g. Ba or Am fluorescence source) vs peak positions to obtain coefficients :math:`p_i`.

Data:

- /detectorCalibration/fit_0001_2_1/fitLine_{ROIi}/horizontal: NXdata plot of peak fit
- /detectorCalibration/fit_0001_2_1/fitLine_{ROIi}/vertical: NXdata plot of peak fit
- /detectorCalibration/fit_0001_2_1/fitParams/fitParamsHD: results of peak fit
        - columns: [height, position, fwhm1, fwhm2, eta, goodness of fit]
        - rows are the peaks
- /detectorCalibration/fit_0001_2_1/fitParams/fitParamsVD: results of peak fit
        - columns: [height, position, fwhm1, fwhm2, eta, goodness of fit]
        - rows are the peaks
- /detectorCalibration/fit_0001_2_1/curveCalibration: NXcollection of energy fit
- /detectorCalibration/fit_0001_2_1/calibCoeffs/calibCoeffsHD: results of energy fit
        - :math:`p_i` stored in `numpy.polyval` order
- /detectorCalibration/fit_0001_2_1/calibCoeffs/calibCoeffsVD: results of energy fit
        - :math:`p_i` stored in `numpy.polyval` order

Task 2: `angleCalibEDD`
+++++++++++++++++++++++

Purpose: calibration of the diffraction angle with known sample

.. math::

    2d\sin\theta = \frac{hc}{E}

Procedure:

- fit peaks in channels
- fit d-spacings vs. peak positions (using energy calibration) to obtain :math:`2\theta`.

Data:

- /angleCalibration/fit_0001_3/fitLine_{ROIi}/horizontal: NXdata plot of peak fit
- /angleCalibration/fit_0001_3/fitLine_{ROIi}/vertical: NXdata plot of peak fit
- /angleCalibration/fit_0001_3/fitParams/fitParamsHD: results of peak fit
        - columns: [height, position, fwhm1, fwhm2, Lorentzian fraction, goodness of fit]
        - rows are the peaks
- /angleCalibration/fit_0001_3/fitParams/fitParamsVD: results of peak fit
        - columns: [height, position, fwhm1, fwhm2, Lorentzian fraction, goodness of fit]
        - rows are the peaks
- /angleCalibration/fit_0001_3/curveAngleCalibration: NXcollection of angle fit
- /angleCalibration/fit_0001_3/calibratedAngle/calibratedAngleHD: results of angle fit
        - :math:`2\theta` of the horizontal detector
- /angleCalibration/fit_0001_3/calibratedAngle/calibratedAngleVD: results of angle fit
        - :math:`2\theta` of the vertical detector

Task 3: `fitEDD`
++++++++++++++++

Purpose: peak positions in channels for unknown sample

Procedure:

- fit peaks in channels

Data:

- /BAIII_AB_4_1_0001_42.1/fit/0000/fitLine_{ROIi}: NXdata plot of fit
- /BAIII_AB_4_1_0001_42.1/fit/0000/fitParams/fitParamsHD: results of peak fit
        - columns: [height, position, fwhm1, fwhm2, eta, goodness of fit]
        - rows are the peaks
- /BAIII_AB_4_1_0001_42.1/fit/0000/fitParams/fitParamsVD: results of peak fit
        - columns: [height, position, fwhm1, fwhm2, eta, goodness of fit]
        - rows are the peaks
- /BAIII_AB_4_1_0001_42.1/positioners:
        - ex, ey, ez, ephi, echi, sy
- /BAIII_AB_4_1_0001_42.1/tthPositionsGroup/peak_{i}: combine positioners and fit results
        - columns: [ex, ey, ez, ephi, echi, sy, delta (azimuth), 2theta, position (in channel), intensity, (FWHM1+FWHM2)/2, Lorentzian fraction, goodness of fit]
        - rows: [horizontal, vertical]

delta is -90 degrees for the horizontal detector (det0) and 0 degrees for the vertical detector (det1).

theta is fixed to zero but should really be the detector angle from `angleCalibEDD`.

Task 4: `coordTransformation`
+++++++++++++++++++++++++++++

Transformation of the coordinates from the motor to the goniometer reference frame.

Data:

- /global/peak_{ROIi}: data in motor frame
        - columns: [ex, ey, ez, ephi, echi, sy, delta (azimuth), 2theta, position (in channel), intensity, (FWHM1+FWHM2)/2, Lorentzian fraction, goodness of fit]
        - rows: [horizontal, vertical]

- /global/inSample_peak_{ROIi}: data in goniometer frame
        - columns: [x, y, z, rot_Z, rot_X, rot_Y, delta, 2theta, position (in channel), intensity, (FWHM1+FWHM2)/2, Lorentzian fraction, goodness of fit]
        - rows: [horizontal, vertical]

x = -ex, y = -ey, z = 9 - ez, rot_Z = ephi, rot_X = echi, rot_Y = sy

Task 5: `regroupPoints`
+++++++++++++++++++++++

Merge data from 2 different sample orientations (2 scans)

Data:

- /coordInSample_Peak_{ROIi}: data in goniometer frame
        - columns: [x, y, z, rot_Z, rot_X, rot_Y, delta, 2theta, position (in channel), intensity, (FWHM1+FWHM2)/2, Lorentzian fraction, goodness of fit]
        - rows: [horizontal scan1, vertical scan1, horizontal scan2, vertical scan2]

Task 6: `preStraind0cstEDD`
+++++++++++++++++++++++++++

Calculate the strain in the measurement direction.

Data:

- /STRAIN_with_d0/peak_{ROIi}/point_00000: data in goniometer frame
        - columns: [x, y, z, rot_minZ, rot_X, rot_minY, delta, 2theta, epsilon_sample, qx, qy, qz]
        - rows: [horizontal scan1, vertical scan1, horizontal scan2, vertical scan2]

epsilon_sample is the strain in the direction of the scattering vector. [qx, qy, qz] is the normalized scattering vector
in cartesian coordinates calculated from spherical coordinates [delta, 2theta] and goniometer rotations rot_minZ, rot_X, rot_minY.

Task 7: `strainStressd0cstEDD`
++++++++++++++++++++++++++++++

Calculate the strain and stress tensors.

Data:

- /peak_{ROIi}/point_00000/position:
        - x, y, z
- /peak_0000/point_00000/strain_tensor_fit:
        - ???
- /peak_0000/point_00000/stress_tensor_fit:
        - ???
