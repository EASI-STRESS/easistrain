
import numpy as np
import h5py
import scipy.optimize

########### Example of the arguments of the main function #############
# fileRead = '/home/esrf/slim/easistrain/easistrain/EDD/Results_ihme10_globstrain_TiC.h5'
# fileSave = '/home/esrf/slim/easistrain/easistrain/EDD/Results_ihme10_strainStressTensor_TiC.h5'
# numberOfPeaks = 5
# XEC = [S1(hkl1), dS2(hkl1), S1(hkl2), dS2(hkl2), S1(hkl3), dS2(hkl3), S1(hkl4), dS2(hkl4), S1(hkl5), dS2(hkl5)]




def deforDirMeas(angles, e11, e22, e33, e23, e13, e12):
    phi = np.radians(angles[:, 0])
    chi = np.radians(angles[:, 1])
    omega = np.radians(angles[:, 2])
    delta = np.radians(angles[:, 3])
    theta = np.radians(0.5*angles[:, 4])
    q1 = (
        (np.cos(theta) * np.cos(chi) * np.sin(delta) * np.sin(phi))
        + (
            np.cos(delta)
            * np.cos(theta)
            * (
                (np.cos(phi) * np.sin(omega))
                - (np.cos(omega) * np.sin(phi) * np.sin(chi))
            )
        )
        - np.sin(theta)
        * ((np.cos(phi) * np.cos(omega)) + (np.sin(phi) * np.sin(chi) * np.sin(omega)))
    )
    q2 = (
        np.cos(delta)
        * np.cos(theta)
        * ((np.cos(phi) * np.cos(omega) * np.sin(chi)) + (np.sin(phi) * np.sin(omega)))
        - (np.cos(theta) * np.cos(phi) * np.cos(chi) * np.sin(delta))
        - (
            np.sin(theta)
            * (
                (np.cos(omega) * np.sin(phi))
                - (np.cos(phi) * np.sin(chi) * np.sin(omega))
            )
        )
    )
    q3 = (
        (np.cos(delta) * np.cos(theta) * np.cos(chi) * np.cos(omega))
        + (np.cos(theta) * np.sin(delta) * np.sin(chi))
        + (np.cos(chi) * np.sin(theta) * np.sin(omega))
    )
    defDirMeas = (
        (e11 * q1 ** 2)
        + (e22 * q2 ** 2)
        + (e33 * q3 ** 2)
        + (2 * e12 * q1 * q2)
        + (2 * e13 * q1 * q3)
        + (2 * e23 * q2 * q3)
    )
    return defDirMeas



def deforDirMeasStress(anglesAndXEC, s11, s22, s33, s23, s13, s12):
    phi = np.radians(anglesAndXEC[:-1, 0])
    chi = np.radians(anglesAndXEC[:-1, 1])
    omega = np.radians(anglesAndXEC[:-1, 2])
    delta = np.radians(anglesAndXEC[:-1, 3])
    theta = np.radians(0.5*anglesAndXEC[:-1, 4])
    S1 = anglesAndXEC[-1,0]
    dS2 = anglesAndXEC[-1,1]
    q1 = (
        (np.cos(theta) * np.cos(chi) * np.sin(delta) * np.sin(phi))
        + (
            np.cos(delta)
            * np.cos(theta)
            * (
                (np.cos(phi) * np.sin(omega))
                - (np.cos(omega) * np.sin(phi) * np.sin(chi))
            )
        )
        - np.sin(theta)
        * ((np.cos(phi) * np.cos(omega)) + (np.sin(phi) * np.sin(chi) * np.sin(omega)))
    )
    q2 = (
        np.cos(delta)
        * np.cos(theta)
        * ((np.cos(phi) * np.cos(omega) * np.sin(chi)) + (np.sin(phi) * np.sin(omega)))
        - (np.cos(theta) * np.cos(phi) * np.cos(chi) * np.sin(delta))
        - (
            np.sin(theta)
            * (
                (np.cos(omega) * np.sin(phi))
                - (np.cos(phi) * np.sin(chi) * np.sin(omega))
            )
        )
    )
    q3 = (
        (np.cos(delta) * np.cos(theta) * np.cos(chi) * np.cos(omega))
        + (np.cos(theta) * np.sin(delta) * np.sin(chi))
        + (np.cos(chi) * np.sin(theta) * np.sin(omega))
    )
    defDirMeasStress = (S1 * (s11 + s22 + s33)) + (dS2 * ((s11 * q1**2)
    + (s22 * q2**2) + (s33 * q3**2) + (2 * s12 * q1 * q2)
    + (2 * s13 * q1 * q3) + (2 * s23 * q2 * q3))
    )
    return defDirMeasStress


def strainStressTensor(
	fileRead,
	fileSave,
	numberOfPeaks,
	XEC
	):
	
	
	h5Save = h5py.File(fileSave, "a")  ## create/append h5 file to save in
	strainTensor = h5Save.create_group('strain_tensor') ## creation of the strain tensor group
	stressTensor = h5Save.create_group('stress_tensor') ## creation of the stress tensor group
	
	with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of the energy calibration of the detectors
		for peakNumber in range(numberOfPeaks):
			peakGroup = strainTensor.create_group(f'peak_{str(peakNumber).zfill(4)}'
			) ## create group for each peak in the strainTensor group
			peakGroupStress = stressTensor.create_group(f'peak_{str(peakNumber).zfill(4)}'
			) ## create group for each peak in the stressTensor group
			for i in range(int(len(list(h5Read[f'STRAIN_with_d0/peak_{str(peakNumber).zfill(4)}']))/2)):
				pointInPeak = h5Read[f'STRAIN_with_d0/peak_{str(peakNumber).zfill(4)}/point_{str(i).zfill(5)}'][()] ##
				upointInPeak = h5Read[f'STRAIN_with_d0/peak_{str(peakNumber).zfill(4)}/uncertainty_point_{str(i).zfill(5)}'][()] ##
				pointInPeakGroup = peakGroup.create_dataset(
				f'point_{str(i).zfill(5)}',
				dtype = 'float64',
				data = np.zeros((1,9),'float64')
				) ## strain tensor in point i
				pointInPeakGroupStress = peakGroupStress.create_dataset(
				f'point_{str(i).zfill(5)}',
				dtype = 'float64',
				data = np.zeros((1,9),'float64')
				) ## stress tensor in point i
				upointInPeakGroup = peakGroup.create_dataset(
				f'uncertainty_point_{str(i).zfill(5)}',
				dtype = 'float64',
				data = np.zeros((1,9),'float64')
				) ## uncertinty on the strain tensor component in point i
				upointInPeakGroupStress = peakGroupStress.create_dataset(
				f'uncertainty_point_{str(i).zfill(5)}',
				dtype = 'float64',
				data = np.zeros((1,9),'float64')
				) ## uncertainty on the stress tensor component in point i
				pointInPeakGroup[0, 0:3] = pointInPeak[0, 0:3] ## coordinate of the point in the sample reference to put in the strain tensor components matrix
				upointInPeakGroup[0, 0:3] = upointInPeak[0, 0:3] ## coordinate of the point in the sample reference to put in the strain tensor components matrix
				pointInPeakGroupStress[0, 0:3] = pointInPeak[0, 0:3] ## coordinate of the point in the sample reference to put in the stress tensor components matrix
				upointInPeakGroupStress[0, 0:3] = upointInPeak[0, 0:3] ## coordinate of the point in the sample reference to put in the stress tensor components matrix
				guessEps11 = np.mean(pointInPeak[pointInPeak[:,9] >= 0.9, 8]) if len(pointInPeak[pointInPeak[:,9] >= 0.9, 8]) > 0 else 0  ## guess of strainxx
				guessEps22 = np.mean(pointInPeak[pointInPeak[:,10] >= 0.9, 8]) if len(pointInPeak[pointInPeak[:,10] >= 0.9, 8]) > 0 else 0 ## guess of strainyy
				guessEps33 = np.mean(pointInPeak[pointInPeak[:,11] >= 0.9, 8]) if len(pointInPeak[pointInPeak[:,11] >= 0.9, 8]) > 0 else 0 ## guess of strainzz
				guessEps23 = np.mean(pointInPeak[pointInPeak[:,10] * pointInPeak[:,11] >= 0.9, 8]) - 0.5*(guessEps22+guessEps33) if len(pointInPeak[pointInPeak[:,10] * pointInPeak[:,11] >= 0.9, 8]) > 0 else 0 ## guess of strainyz
				guessEps13 = np.mean(pointInPeak[pointInPeak[:,9] * pointInPeak[:,11] >= 0.9, 8]) - 0.5*(guessEps11+guessEps33) if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,11] >= 0.9, 8]) > 0 else 0 ## guess of strainxz
				guessEps12 = np.mean(pointInPeak[pointInPeak[:,9] * pointInPeak[:,10] >= 0.9, 8]) - 0.5*(guessEps11+guessEps22) if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,10] >= 0.9, 8]) > 0 else 0 ## guess of strainxz
				guessSig11 = (1/XEC[2 * peakNumber + 1]) * (guessEps11
				+ ((- XEC[2 * peakNumber] / (XEC[2 * peakNumber + 1] + 3 * XEC[2 * peakNumber])) * (guessEps11 + guessEps22 + guessEps33))) ## guess of sigmaxx
				guessSig22 = (1/XEC[2 * peakNumber + 1]) * (guessEps22
				+ ((- XEC[2 * peakNumber] / (XEC[2 * peakNumber + 1] + 3 * XEC[2 * peakNumber])) * (guessEps11 + guessEps22 + guessEps33))) ## guess of sigmayy
				guessSig33 = (1/XEC[2 * peakNumber + 1]) * (guessEps33
				+ ((- XEC[2 * peakNumber] / (XEC[2 * peakNumber + 1] + 3 * XEC[2 * peakNumber])) * (guessEps11 + guessEps22 + guessEps33))) ## guess of sigmazz
				guessSig23 = (1/XEC[2 * peakNumber + 1]) * guessEps23 ## guess of sigmayz
				guessSig13 = (1/XEC[2 * peakNumber + 1]) * guessEps13 ## guess of sigmaxz
				guessSig12 = (1/XEC[2 * peakNumber + 1]) * guessEps12 ## guess of sigmaxy
				strainTensorInitialGuess = [guessEps11, guessEps22, guessEps33,
				guessEps23, guessEps13, guessEps12] ## initial guess of the strain tensor
				stressTensorInitialGuess = [guessSig11, guessSig22, guessSig33,
				guessSig23, guessSig13, guessSig12] ## initial guess of the stress tensor
				boundsEps11 = [-np.inf if len(pointInPeak[pointInPeak[:,9] >= 0.1, 8]) > 0 else -10**-10, 
				np.inf if len(pointInPeak[pointInPeak[:,9] >= 0.1, 8]) > 0 else 10**-10] ## bounds on epsilon11
				boundsEps22 = [-np.inf if len(pointInPeak[pointInPeak[:,10] >= 0.1, 8]) > 0 else -10**-10, 
				np.inf if len(pointInPeak[pointInPeak[:,10] >= 0.1, 8]) > 0 else 10**-10] ## bounds on epsilon22
				boundsEps33 = [-np.inf if len(pointInPeak[pointInPeak[:,11] >= 0.1, 8]) > 0 else -10**-10,
				np.inf if len(pointInPeak[pointInPeak[:,11] >= 0.1, 8]) > 0 else 10**-10] ## bounds on epsilon33
				boundsEps23 = [-np.inf if len(pointInPeak[pointInPeak[:,10] * pointInPeak[:,11] >= 0.1, 8]) > 0 else -10**-10,
				np.inf if len(pointInPeak[pointInPeak[:,10] * pointInPeak[:,11] >= 0.1, 8]) > 0 else 10**-10] ## bounds on epsilon23
				boundsEps13 = [-np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,11] >= 0.1, 8]) > 0 else -10**-10,
				np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,11] >= 0.1, 8]) > 0 else 10**-10] ## bounds on epsilon13
				boundsEps12 = [-np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,10] >= 0.1, 8]) > 0 else -10**-10,
				np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,10] >= 0.1, 8]) > 0 else 10**-10] ## bounds on epsilon12
				boundsSig11 = [-np.inf, np.inf] ## bounds on sigma11
				boundsSig22 = [-np.inf, np.inf]  ## bounds on sigma22
				boundsSig33 = [-np.inf , np.inf] ## bounds on sigma33
				boundsSig23 = [-np.inf if len(pointInPeak[pointInPeak[:,10] * pointInPeak[:,11] >= 0.1, 8]) > 0 else -10**-10,
				np.inf if len(pointInPeak[pointInPeak[:,10] * pointInPeak[:,11] >= 0.1, 8]) > 0 else 10**-10] ## bounds on sigma23
				boundsSig13 = [-np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,11] >= 0.1, 8]) > 0 else -10**-10,
				np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,11] >= 0.1, 8]) > 0 else 10**-10] ## bounds on sigma13
				boundsSig12 = [-np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,10] >= 0.1, 8]) > 0 else -10**-10,
				np.inf if len(pointInPeak[pointInPeak[:,9] * pointInPeak[:,10] >= 0.1, 8]) > 0 else 10**-10] ## bounds on sigma12
				#print(strainTensorInitialGuess)
				#print(stressTensorInitialGuess)
				#print(([boundsEps11[0],boundsEps22[0],boundsEps33[0],boundsEps23[0],boundsEps13[0],boundsEps12[0]],
				#[boundsEps11[1],boundsEps22[1],boundsEps33[1],boundsEps23[1],boundsEps13[1],boundsEps12[1]]))
				#print(([boundsSig11[0],boundsSig22[0],boundsSig33[0],boundsSig23[0],boundsSig13[0],boundsSig12[0]],
				#[boundsSig11[1],boundsSig22[1],boundsSig33[1],boundsSig23[1],boundsSig13[1],boundsSig12[1]]))
				strainTensorComponents, covarStrains = scipy.optimize.curve_fit(
				f = deforDirMeas,
				xdata = pointInPeak[:, 3:8],
				ydata = pointInPeak[:, 8],
				p0 = strainTensorInitialGuess,
				sigma = upointInPeak[:, 8],
				bounds = ([boundsEps11[0],boundsEps22[0],boundsEps33[0],boundsEps23[0],boundsEps13[0],boundsEps12[0]],
				[boundsEps11[1],boundsEps22[1],boundsEps33[1],boundsEps23[1],boundsEps13[1],boundsEps12[1]])
				) ## fit of the strain tensor
				stressTensorComponents, covarStress = scipy.optimize.curve_fit(
				f = deforDirMeasStress,
				xdata = np.append(pointInPeak[:, 3:8],np.array([[XEC[2 * peakNumber],XEC[2 * peakNumber + 1],0,0,0]]), axis = 0),
				ydata = pointInPeak[:, 8],
				p0 = stressTensorInitialGuess,
				sigma = upointInPeak[:, 8],
				bounds = ([boundsSig11[0],boundsSig22[0],boundsSig33[0],boundsSig23[0],boundsSig13[0],boundsSig12[0]],
				[boundsSig11[1],boundsSig22[1],boundsSig33[1],boundsSig23[1],boundsSig13[1],boundsSig12[1]])
				) ## fit of the stress tensor
				pointInPeakGroup[0, 3:9] = strainTensorComponents  ## strain tensor components
				upointInPeakGroup[0, 3:9] = np.sqrt(np.diag(covarStrains)) if np.sum(covarStrains) != np.inf else np.mean(upointInPeak[:, 8]) ## uncertainty on the strain tensor components
				pointInPeakGroupStress[0, 3:9] = stressTensorComponents  ## stress tensor components
				upointInPeakGroupStress[0, 3:9] = np.sqrt(np.diag(covarStress)) if np.sum(covarStress) != np.inf else (1/(XEC[2 * peakNumber + 1] + XEC[2 * peakNumber])) * np.mean(upointInPeak[:, 8]) ## uncertainty on the strain tensor components
				
	h5Save.close()
	return
