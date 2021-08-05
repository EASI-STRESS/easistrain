
from func_lattice_d_spacing_param import *
import numpy as np
import sys
import time

time0 = time.time()

if (len(sys.argv) != 7 or sys.argv[1] == 'help'):
	print('#### Arguments to be given: 1)root // 2)h5file // 3)path of poni file // 4)h // 5)k // 6)l')
else:
	print('\n\n#*#*#*#*#*#*#*#*#* Calculation of lattice and d-spacing parameters procedure started #*#*#*#*#*#*#*#*#*#*#*\n\n')
	lattice_param(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
	print('\n\n#*#*#*#*#*#*#*#*#* Calculation of lattice and d-spacing parameters procedure finished #*#*#*#*#*#*#*#*#*#*#*\n\n')

time1 = time.time()

print('total time: ' + str(time1-time0) + ' seconds')
