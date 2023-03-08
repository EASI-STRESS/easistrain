import os
import matplotlib.pyplot as plt
import h5py


def V1D(root_data, h5file, dataset):
    """Example of the arguments of the main function"""
    filename = os.path.join(root_data, "Result", h5file)
    with h5py.File(filename, "r") as f:
        path = f[dataset][()]
    plt.figure(figsize=(10, 8))
    plt.plot(path[:, 0], path[:, 1:])
    plt.show()


def interV1D(root_data, h5file, dsetX, dsetY, Xlabel, Ylabel):
    with h5py.File(os.path.join(root_data, h5file), "r") as f:
        path1 = f[dsetX][()]
        path2 = f[dsetY][()]
    plt.figure(figsize=(10, 8))
    plt.plot(path1, path2)
    plt.xlabel(Xlabel, family="sans-serif", fontsize=28)
    plt.ylabel(Ylabel, family="sans-serif", fontsize=28)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    plt.show()
