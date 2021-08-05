

import h5py
import numpy as np


### This function get the image matrix from the h5 file and convert it to float64 ###
### root_data: the path of the folder where the h5 file is saved
### h5file: The name of the h5 file from which the matrix of the image will be extracted
### scan: The name of the group on which the concerned measurement are saved
### detector name: The name of the detector

def get_image_matrix(root_data, h5file, scan, detector_name):
	r_h5file = h5py.File(root_data + '/' + h5file,'r')
	image = r_h5file['/' + scan + '/measurement/' + detector_name]
	if (np.ndim(image) == 2):
		print(np.shape(image))
		image_matrix = np.float64(image)
	else:
		print(np.shape(image))
		print('### The image matrix is not a 2D squared matrix')
		image_matrix = np.float64(image[0,:,:])
	return image_matrix