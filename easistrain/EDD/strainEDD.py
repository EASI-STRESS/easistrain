
import h5py
import numpy as np


#fileRead = '/home/esrf/slim/easistrain/easistrain/EDD/Results_ihme10_test_TiC.h5'
#fileSave = '/home/esrf/slim/easistrain/easistrain/EDD/Results_ihme10_test_TiC.h5'
#numberOfPeaks = 5


def strainEDD(
	fileRead,
	fileSave,
	numberOfPeaks
	
	):


	
	with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
		scanList = list(h5Read.keys())  ## list of the scans
		lengthCounter = 0	
		for scan in scanList:
			shapeDsetPeak = h5Read[f'{scan}/tthPositionsGroup/peak_0000'].shape[0]
			lengthCounter = lengthCounter + shapeDsetPeak ## shape of the matrix on which all the points for one peak will be saved
		
	h5Save = h5py.File(fileSave, "a")  ## create/append h5 file to save in
	globalGroup = h5Save.create_group(
    "global",
    )  ## Creation of the global group in which all peaks positions of all the points will be put
	for peakNumber in range(numberOfPeaks):
		globalPeak = globalGroup.create_dataset(
		f"peak_{str(peakNumber).zfill(4)}",
		dtype="float64",
		data = np.zeros((lengthCounter,13),'float64')
		) ## creation of the dataset for each peak
		arrowsCounter = 0
		for scan in scanList:
			with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
				shapeDsetPeak = h5Read[f'{scan}/tthPositionsGroup/peak_0000'].shape[0]
				globalPeak[arrowsCounter:arrowsCounter + shapeDsetPeak,:] = h5Read[f'{scan}/tthPositionsGroup/peak_{str(peakNumber).zfill(4)}'][()]
				arrowsCounter = arrowsCounter + shapeDsetPeak
	for peakNumber in range(numberOfPeaks):
		peakGroup = globalGroup.create_group(
		f'pointsPerPeak_{str(peakNumber).zfill(4)}'
		) ## create a group per peak in the global group
		pointsCounter = 0
		matToFilter = globalGroup[f'peak_{str(peakNumber).zfill(4)}'][()]
		matCheck = np.array([[]])
		for i in range(len(matToFilter)):
			matFirstFilter = matToFilter[matToFilter[:,0] == matToFilter[i,0],
			:] ## filtering based on x coordinate of the measurement point
			matSecFilter = matFirstFilter[matFirstFilter[:,1] == matFirstFilter[0,1], 
			:] ## filtering based on y coordinate of the measurement point
			matThirdFilter = matSecFilter[matSecFilter[:,2] == matSecFilter[0,2], 
			:] ## filtering based on z coordinate of the measurement point
			if pointsCounter == 0:
				peakGroup.create_dataset(
				f'point_{str(pointsCounter).zfill(5)}',
				dtype = 'float64',
				data = matThirdFilter
				)
				pointsCounter = pointsCounter + 1
				matCheck = np.append(matCheck, matThirdFilter[0,:3])
				matCheck = np.reshape(matCheck, (int(len(matCheck)/3),3))
			if 	pointsCounter > 0:
				check = np.array(())
				for kk in range(len(matCheck)):
					check = np.append(check,np.sum(matThirdFilter[0,:3] == matCheck[kk,:3]))
				if np.sum(check[:] == 3) == 0:
					peakGroup.create_dataset(
					f'point_{str(pointsCounter).zfill(5)}',
					dtype = 'float64',
					data = matThirdFilter
					)
					pointsCounter = pointsCounter + 1
					matCheck = np.append(matCheck, matThirdFilter[0,:3])
					matCheck = np.reshape(matCheck, (int(len(matCheck)/3),3))
				else :
					pointsCounter = pointsCounter
	## Do the same with uncertainty as i did wit peaks and verify for the peaks
	## modifications				
					
					
			
	
	
	
	h5Save.close()		
	return


#peakInfo = h5Read[
#f'{scan}/tthPositionsGroup/peak_{str(peak).zfill(4)}'][()] ## info of the peak (I,pos,FWHM,shape factor), coordinates of the point and angles
#uncertaintyPeakInfo = h5Read[
#f'{scan}/tthPositionsGroup/uncertaintyPeak_{str(peak).zfill(4)}'][()] ## uncertainty
