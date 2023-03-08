import os
import h5py
import numpy


def get_image_matrix(root_data, h5file, scan, detector_name):
    """Get the image matrix from the h5 file and convert it to float64.

    :param root_data: the path of the folder where the h5 file is saved
    :param h5file: The name of the h5 file from which the matrix of the image will be extracted
    :param scan: The name of the group on which the concerned measurement are saved
    :param detector name: The name of the detector
    """
    with h5py.File(os.path.join(root_data, h5file), "r") as f:
        image = f[f"/{scan}/measurement/{detector_name}"]
        if numpy.ndim(image) == 2:
            print(numpy.shape(image))
            image_matrix = numpy.float64(image)
        else:
            print(numpy.shape(image))
            print("### The image matrix is not a 2D squared matrix")
            image_matrix = numpy.float64(image[0, :, :])
        return image_matrix
