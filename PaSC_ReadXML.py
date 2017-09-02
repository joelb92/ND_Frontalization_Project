import os
#import numpy as np
#import matplotlib.pyplot as plt
import sys
from shutil import copyfile

f1 = open('pasc_video_control.xml','r')
f2 = open('pasc_control_names.txt','w')
cnt = 0
for line in f1:
	if 'file-name' in line:
		fl = line.split('file-name=')[1].split('"')[1].split('.')[0]
		f2.write(fl+'\n')
		cnt += 1
print cnt
f1.close()
f2.close()


