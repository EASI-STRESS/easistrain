
from func_integration_2D import *
import numpy as np
import sys
import time

time0 = time.time()

if (len(sys.argv) != 11 or sys.argv[1] == 'help'):
	print('#### Arguments to be given: 1)root // 2)h5file // 3)scan // 4)detector_name // 5)path of poni file // 6)npt_rad // 7)npt_azim // 8)unit(2th_deg or 2th_rad or q_A^-1...) //9)dark image if exist if not give 0//10)mask image if exist if not give 0')
else:
	integration_2D(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])

time1 = time.time()

print('total time: ' + str(time1-time0) + ' seconds')
