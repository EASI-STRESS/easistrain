import os
import numpy
import fabio
from .func_save_edf_image import save_edf_image


def mask(root, image, int_min, int_max, root_save, mask_name, extension):
    """Generate automatically a mask for an image. It masks the pixels lower than int_min and bigger than int_max

    :param root: the path of the folder of the image for which the mask will be generated
    :param image: The name of the image for which the mask will be generated
    :param int_min: minimum intensity below which the pixels have to be masked
    :param int_max: maximum intensity above which the pixels have to be masked
    :param root_save: The path of the folder pn which the mask have to be saved
    :param mask_name: The name of the mask
    :param extension: The format of the mask image (edf, tif, ...)"""
    image_matrix = fabio.open(os.path.join(root, image + "." + extension)).data
    mask_matrix = numpy.zeros_like(image_matrix)
    for i in range(numpy.shape(image_matrix)[0]):
        for j in range(numpy.shape(image_matrix)[1]):
            if image_matrix[i, j] < numpy.float64(int_min) or image_matrix[
                i, j
            ] > numpy.float64(int_max):
                mask_matrix[i, j] = 1
    save_edf_image(root_save, mask_name, extension, mask_matrix)
