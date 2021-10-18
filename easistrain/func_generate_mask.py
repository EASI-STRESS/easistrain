import fabio
import numpy as np
from fabio.edfimage import edfimage
from func_save_edf_image import *

### This function generate automatically a mask for an image. It masks the pixels lower than int_min and bigger than int_max
### root: the path of the folder of the image for which the mask will be generated
### image: The name of the image for which the mask will be generated
### int_min: minimum intensity below which the pixels have to be masked
### int_max: maximum intensity above which the pixels have to be masked
### root_save: The path of the folder pn which the mask have to be saved
### mask_name: The name of the mask
### extension: The format of the mask image (edf, tif, ...)


def mask(root, image, int_min, int_max, root_save, mask_name, extension):
    image_matrix = fabio.open(root + "/" + image + "." + extension).data
    mask_matrix = np.zeros_like(image_matrix)
    for i in range(np.shape(image_matrix)[0]):
        for j in range(np.shape(image_matrix)[1]):
            if image_matrix[i, j] < np.float64(int_min) or image_matrix[
                i, j
            ] > np.float64(int_max):
                mask_matrix[i, j] = 1
    save_edf_image(root_save, mask_name, extension, mask_matrix)
    return
