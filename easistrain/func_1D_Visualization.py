
import matplotlib.pyplot as plt
import numpy as np
import h5py


def V1D(root_data, h5file, dataset):
	r_h5file = h5py.File(root_data + '/' + 'Results' + '_' + h5file,'r') ### Read the h5 file
	path = r_h5file[dataset]
	plt.figure(figsize=(10,8))
	plt.plot(path[:,0],path[:,1:])
	plt.show()
	return

