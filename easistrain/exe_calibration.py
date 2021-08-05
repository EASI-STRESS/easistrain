
from func_get_image_matrix import * 
from func_save_edf_image import *
from func_generate_mask import *
import numpy as np
import sys

### This routine appy three function:
### 1) get_image_matrix: gets the matrix of the calibration image from h5 file
### 2) save_edf_image: saves the calibration
### 3) mask: generates a mask for the image and saves it in the sae folder of the calibration image 

if len(sys.argv) != 10:  ### Verification of the number of introduced arguments in the command line
	print('####  The number of arguments is not correct')
	print('#### Arguments to be given: 1)root_data // 2)h5file // 3)scan // 4)detector_name // 5)root_save // 6)save_image_name // 7)extension // 8)int_min // 9)int_max')
	quit()
else:
	print('#### Enjoy the calibration')
	image_matrix = get_image_matrix(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])  ### Usage of the get_image_matrix function (see function)
	print('#### image matrix OK')
	print('dimension of the matirx image =' + str(np.shape(image_matrix)))
	save_edf_image(sys.argv[5], sys.argv[6], sys.argv[7], image_matrix) ### Usage of the function save_edf_image to save the calibration image
	print('#### The calibrant image was saved in edf format in:' + ' ' + sys.argv[5] + '/' + sys.argv[6] + '.' + sys.argv[7])
	mask(sys.argv[5], sys.argv[6], sys.argv[8], sys.argv[9], sys.argv[5], sys.argv[6] + '_mask', sys.argv[7]) ### Usage of the mask image to generate a mask for the calibration image
	print('#### The mask for the calibrant image was saved in:' + ' ' + sys.argv[5] + '/' + sys.argv[6] + '_mask' + '.' + sys.argv[7])

print('#### Now we will quit python and let you do the calibration')
print('#### When python is quitted please tape:')
print('#### pyFAI-calib -e (energy in keV)  -p (pixel size in micron) -P 0.95 (Polarization) --npt (nb of points in 1D integrated pattern) --npt-azim (nb of azimuthal sectors) -c (calibrant) -m (mask image) calib_image')
