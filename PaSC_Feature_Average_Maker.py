import numpy as np
import os
import sys
#from scipy import spatial

fold = sys.argv[1]
f1 = open(fold+'/pasc_control_feat.txt','r')
f2 = open(fold+'/pasc_control_feat_Averaged.txt','w')

feat_len = 4096  # length of feature vector (4096 for VGG fc7)
prev_cls = 0
cnt = 0
for line in f1:
	fl = line.split(',')
	fll = fl[0].split('-')
	cls = fll[0]
	if prev_cls == 0:
		cls_cnt = 1
		prev_cls = cls
		temp_lst = [0.0]*feat_len
		for i in range(1,len(fl)):
			temp_lst[i-1] += (float(fl[i]))
	elif prev_cls == cls:
		cls_cnt += 1
		for i in range(1,len(fl)):
			temp_lst[i-1] += (float(fl[i]))
	elif prev_cls != cls:
		temp_lst[:] = [x / cls_cnt for x in temp_lst]
		f2.write(prev_cls)
		for j in range(0,len(temp_lst)):
			f2.write(',')
			f2.write(str(temp_lst[j]))
		f2.write('\n')
		cls_cnt = 1
		prev_cls = cls
		temp_lst = [0.0]*feat_len
		for i in range(1,len(fl)):
			temp_lst[i-1] += (float(fl[i]))
	cnt += 1
	print cnt

temp_lst[:] = [x / cls_cnt for x in temp_lst]
f2.write(prev_cls)
for j in range(0,len(temp_lst)):
	f2.write(',')
	f2.write(str(temp_lst[j]))
f2.write('\n')

f1.close()
f2.close()

