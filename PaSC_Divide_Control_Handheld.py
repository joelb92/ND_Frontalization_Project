import os
#import numpy as np
#import matplotlib.pyplot as plt
import sys
from distutils.dir_util import copy_tree

str1 = '/scratch365/sbanerj1/PaSC_dlib/handheld'
str2 = '/scratch365/sbanerj1/PaSC_Old_Vito_Asym/asym'
str3 = '/scratch365/sbanerj1/PaSC_Old_Vito_Asym/handheld'
#str4 = '/scratch365/sbanerj1/PaSC_ZR/PaSC_zr_hassner/asym/handheld'

d_list = []
cnt = 0
for dirpath, dirnames, filenames in os.walk(str1):
	for d in dirnames:
		d_list.append(d)
		cnt += 1
print cnt

cnt = 0

for dirpath, dirnames, filenames in os.walk(str2):
	for d in dirnames:
		if d in d_list:
			os.makedirs(str3+'/'+d)
			copy_tree(str2+'/'+d, str3+'/'+d)
			cnt += 1
print cnt

