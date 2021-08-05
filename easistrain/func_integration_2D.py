
import numpy as np
import pyFAI
import h5py
import fabio

### This function integrates a 2D image using integrate2D pyfai's function and save the results in a h5file named Results_name of the h5 file containing the image
### 1) It looks on the image on th h5 file
### 2) creates a mask based on the int_max and int_min the user gives
### 3) Integrate the 2D image in number of given sectors (npt_azim) and with a given number of points in the radial direction (npt_rad)
### 4) saves the results in a h5file named Results_name of the h5 file containing the image
### Definition of the function inputs:
### root_data: the path of the file where the h5 file is saved
### h5file: the name of the h5 file where the images are saved
### scan: The name of the group (in h5file) on which the concerned measurement are saved 
### detector name: name of rhe used detector
### poni_file: the path of the poni file which define the experimental geometry (ai)
### int_min: minimum intensity below which the pixels have to be masked (deleted)
### int_max: maximum intensity above which the pixels have to be masked (deleted)
### npt_rad: number of points in the radial direction (number of 2theta points)
### npt_azim: number of sectors for the integration in the azimutha direction (around the ring)
### x_unit: the unity of the x axis (2th_deg or 2th_rad or q_A^-1...)
### im_dark: the name of the dark image (if no dark image exist please give 0 as argument)
### im_mask: the name of the mask image (if no mask image exist please give 0 as argument)

def integration_2D(root_data, h5file, scan, detector_name, poni_file, npt_rad, npt_azim, x_unit, im_dark, im_mask):
	print(im_dark)
	fh5_save = h5py.File(root_data + '/' + 'Results' + '_' + h5file,'a') ### Create th file in which will be saved the results (integration, ...)
	level_1 = fh5_save.create_group(scan) ### Create the group on which will be putted the data of integration
	level_1_subg_1 = level_1.create_group('raw_integration_2D')
	level_1_subg_3 = level_1_subg_1.create_group('Integration_parameter')
	level_1_subg_3.create_dataset('npt_azim', dtype ='f', data = int(npt_azim))
	level_1_subg_3.create_dataset('npt_rad', dtype ='f', data = int(npt_rad))
	ai = pyFAI.load(poni_file)  ### Load the poni file describing the integration geometry
	r_h5file = h5py.File(root_data + '/' + h5file,'r') ### Read the h5 file
	image_nc = r_h5file['/' + scan + '/measurement/' + detector_name] ### Read the image or images
	image = np.float64(image_nc) # convert the image matrix from int32 to float
	print('### The matrix image is a ' +  str(np.ndim(image)) + 'D matrix')
	if (np.ndim(image) == 2): # do a test on the dimension of the image matrix
		rslt_matrix_cts = np.zeros((int(npt_rad), int(npt_azim) + 1), float)
		rslt_matrix_chi = np.zeros((int(npt_azim)), float)
		rslt_matrix_tth = np.zeros((int(npt_rad)), float)
		if im_mask == '0':
		    mask_mat = np.zeros((np.shape(image)[0], np.shape(image)[1]), float) # creating a zeroes matrix as no mask will be used in this case
		    print('### No mask was used for the integration')
		else:
		    mask_mat = fabio.open(im_mask).data
		    print('### The image: ' + im_mask + ' ' + 'was used as mask')
		#print('### Do not worry I am creating the mask')
		# for ix in range(np.shape(image)[0]): # loop on the image dimension and generation of the mask
		# 	for iy in range(np.shape(image)[1]):
		# 		if (image[ix,iy] < np.float(int_min) or image[ix,iy] > np.float(int_max)):
		# 			mask_mat[ix,iy] = 1 
		# 			#print(str(ix) + '/' + str(iy))
		#print('### Hola, I finished creating the mask')
		print('### Integration started with a'+ ' ' + npt_rad + ' ' + 'steps in tth and with a' + ' ' + npt_azim + ' ' + 'sectors')
		if im_dark == '0': # integration of the image
			cts, tth, chi = ai.integrate2d(image, int(npt_rad), int(npt_azim), correctSolidAngle = True, polarization_factor=0.95, unit=x_unit, method="splitpixel", mask=mask_mat)
		else:
			cts, tth, chi = ai.integrate2d(image, int(npt_rad), int(npt_azim), correctSolidAngle = True, polarization_factor=0.95, unit=x_unit, method="splitpixel", mask=mask_mat, dark = im_dark)
		print('### Integration copmleted')
		#print(np.shape(tth))
		#print(np.shape(cts))
		#print(np.shape(rslt_matrix_cts))
		rslt_matrix_cts[:,0] = tth[:]
		rslt_matrix_chi[:] = chi[:]
		rslt_matrix_tth[:] = tth[:]
		print('### Saving in h5file started')
		for ichi in range(np.shape(chi)[0]):
			rslt_matrix_cts[:,ichi+1] = cts[ichi,:]
		level_1_subg_2 = level_1_subg_1.create_group('image_00000')
		level_1_subg_2.create_dataset('tth_vs_cts', dtype ='f', data = rslt_matrix_cts)
		level_1_subg_2.create_dataset('chi', dtype ='f', data = rslt_matrix_chi)
		level_1_subg_2.create_dataset('tth', dtype ='f', data = rslt_matrix_tth)
		print('### Saving in h5file completed')
		print('### Hope to see you again')
	else:
		#print('### The matrix image is a ' +  str(np.ndim(image)) + 'D matrix')
		if im_mask == '0':
		    mask_mat = np.zeros((np.shape(image)[0], np.shape(image)[1]), float) # creating a zeroes matrix to be filled later with the mask
		    print('### No mask was used for the integration')
		else:
		    mask_mat = fabio.open(im_mask).data
		    print('### The image: ' + im_mask + ' ' + 'was used as mask')
		for i in range(np.shape(image)[0]):
			print('*#*#*#*#* Processing of image_' + '%.0f'%i + '*#*#*#*#*')
			rslt_matrix_cts = np.zeros((int(npt_rad), int(npt_azim) + 1), float)
			rslt_matrix_chi = np.zeros((int(npt_azim)), float)
			rslt_matrix_tth = np.zeros((int(npt_rad)), float)
			#mask_mat = np.zeros((np.shape(image)[0], np.shape(image)[1]), float)
			#print('### Do not worry I am creating the mask')
			# for ix in range(np.shape(image)[1]):
			# 	for iy in range(np.shape(image)[2]):
			# 		if (image[i,ix,iy] < np.float(int_min) or image[i,ix,iy] > np.float(int_max)):
			# 			mask_mat[ix,iy] = 1
			#print('### Hola, I finished creating the mask')
			print('### Integration started with a'+ ' ' + npt_rad + ' ' + 'steps in tth and with a' + ' ' + npt_azim + ' ' + 'sectors' )
			if im_dark == '0':
				cts, tth, chi = ai.integrate2d(image[i], int(npt_rad), int(npt_azim), correctSolidAngle = True, polarization_factor=0.95, unit = x_unit, method = 'splitpixel', mask = mask_mat)
			else:	
				cts, tth, chi = ai.integrate2d(image[i], int(npt_rad), int(npt_azim), correctSolidAngle = True, polarization_factor=0.95, unit = x_unit, method = 'splitpixel', mask = mask_mat, dark = im_dark)
			print('### Integration copmleted')
			rslt_matrix_cts[:,0] = tth[:]
			rslt_matrix_chi[:] = chi[:]
			rslt_matrix_tth[:] = tth[:]
			print('### Saving in h5file started')
			for ichi in range(np.shape(chi)[0]):
				rslt_matrix_cts[:,ichi+1] = cts[ichi,:]
			if (i < 10):
				level_1_subg_2 = level_1_subg_1.create_group('image' + '_0000' + '%.0f'%i)
			if ((i>9) and (i<100)):
				level_1_subg_2 = level_1_subg_1.create_group('image' + '_000' + '%.0f'%i)
			if ((i>99) and (i<1000)):
				level_1_subg_2 = level_1_subg_1.create_group('image' + '_00' + '%.0f'%i)
			if ((i>999) and (i<10000)):
				level_1_subg_2 = level_1_subg_1.create_group('image' + '_0' + '%.0f'%i)
			if ((i>9999) and (i<100000)):
				level_1_subg_2 = level_1_subg_1.create_group('image' + '_' + '%.0f'%i)
			level_1_subg_2.create_dataset('tth_vs_cts', dtype ='f', data = rslt_matrix_cts)
			level_1_subg_2.create_dataset('chi', dtype ='f', data = rslt_matrix_chi)
			level_1_subg_2.create_dataset('tth', dtype ='f', data = rslt_matrix_tth)
			print('### Saving in h5file completed')
		return