import matplotlib.pyplot as plt
import numpy as np
import h5py


def V1D(root_data, h5file, dataset):
    r_h5file = h5py.File(
        root_data + "/" + "Results" + "_" + h5file, "r"
    )  ### Read the h5 file
    path = r_h5file[dataset]
    plt.figure(figsize=(10, 8))
    plt.plot(path[:, 0], path[:, 1:])
    plt.show()
    return


def interV1D(root_data, h5file, dsetX, dsetY, Xlabel, Ylabel):
    r_h5file = h5py.File(root_data + "/" + h5file, "r")  ### Read the h5 file
    path1 = r_h5file[dsetX]
    path2 = r_h5file[dsetY]
    plt.figure(figsize=(10, 8))
    plt.plot(path1[:], path2[:])
    plt.xlabel(Xlabel, family="sans-serif", fontsize=28)
    plt.ylabel(Ylabel, family="sans-serif", fontsize=28)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    plt.show()
    return
