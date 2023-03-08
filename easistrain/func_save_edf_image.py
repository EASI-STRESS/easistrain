import os
from fabio.edfimage import edfimage


def save_edf_image(root_save, save_img_name, extension, img_matrix):
    """
    :param root_save: path of the folder where the image will be saved
    :param save_img_name: name of the image
    :param extension: file extension
    :param img_matrix: the matrix to be saved as image
    """
    iheader = {}
    edfimg = edfimage(data=img_matrix, header=iheader)
    filename = os.path.join(root_save, f"{save_img_name}.{extension}")
    edfimg.write(filename)
