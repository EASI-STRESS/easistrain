
from func_fitting_peaks import *
import numpy as np
import sys
import time

time0 = time.time()

if (len(sys.argv) != 11 or sys.argv[1] == 'help'):
	print('#### Arguments to be given: 1)root // 2)h5file // 3) scan // 4)tth_min // 5)tth_max // 6)bgd_left // 7)bgd_right // 8)fitting function:PVII, gauss, Lorentz or Ps // 9)hkl // 10)Rp:godness of fit')
else:
	print('\n\n#*#*#*#*#*#*#*#*#* Fitting procedure started #*#*#*#*#*#*#*#*#*#*#*\n\n')
	fit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
	print('\n\n#*#*#*#*#*#*#*#*#* Fitting procedure completed #*#*#*#*#*#*#*#*#*#*#*\n\n')
	print('\n\n#*#*#*#*#*#*#*#*#* Cleaning procedure started #*#*#*#*#*#*#*#*#*#*#*\n\n')
	clean_fit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
	print('\n\n#*#*#*#*#*#*#*#*#* Cleaning procedure completed #*#*#*#*#*#*#*#*#*#*#*\n\n')
	

time1 = time.time()

print('total time: ' + str(time1 - time0) + ' seconds')