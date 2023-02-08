import numpy as np
import fabio
import configparser
from fabio.edfimage import edfimage
from fabio.tifimage import tifimage


### This function save the calibrant image as edf extension ###
### root_save: path of the folder where the image will be saved
### save_img_name: name of the image
### extension: format (extension) of the image (edf, tif, ...)
### img_matrix: the matrix to be saved as image
def save_edf_image(root_save, save_img_name, extension, img_matrix):
    iheader = {}
    edfimg = edfimage(data=img_matrix, header=iheader)
    edfimg.write(root_save + "/" + save_img_name + "." + extension)
    return
