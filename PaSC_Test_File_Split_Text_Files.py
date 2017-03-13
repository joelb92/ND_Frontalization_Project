import os
#import numpy as np
#import matplotlib.pyplot as plt
import sys

str1 = '/scratch365/sbanerj1/PaSC_Old_Vito_Asym/handheld'
c = 0
cnt_img = 1
cnt_file = 1
f1 = open('/scratch365/sbanerj1/PaSC_Old_Vito_Asym/pasc_handheld'+str(cnt_file)+'.txt','w')
for dirpath, dirnames, filenames in os.walk(str1):
	for file1 in [f for f in filenames if f.endswith(".jpg")]:
		if cnt_img <= 10000:
			f1.write(file1+'\n')
			cnt_img += 1
			c += 1
		else:
			f1.close()
			cnt_file += 1
			cnt_img = 1
			f1 = open('//scratch365/sbanerj1/PaSC_Old_Vito_Asym/pasc_handheld'+str(cnt_file)+'.txt','w')
			f1.write(file1+'\n')
			cnt_img += 1
			c += 1
print c
f1.close()
	
'''
str1 = '/scratch365/sbanerj1/PaSC_2DAligned/test'
str2 = '/scratch365/sbanerj1/PaSC_2DAligned/control_testimages_nobaddetections.csv'

f2 = open(str2,'r')

det_list = []
c = 0
for line in f2:
	#fl = line.split("/")
	fl = line.split("'")
	det_list.append(fl[1])
	#print fl[1]
	c += 1
print c
print det_list[0], det_list[1]

#f1 = open('rand.txt','w')

c = 0
cnt_img = 1
cnt_file = 1
f1 = open('/scratch365/sbanerj1/PaSC_2DAligned/pasc_control'+str(cnt_file)+'.txt','w')
for dirpath, dirnames, filenames in os.walk(str1):
	for file1 in [f for f in filenames if f.endswith(".jpg")]:
		if file1 in det_list:
			if cnt_img <= 10000:
				f1.write(file1+'\n')
				cnt_img += 1
				c += 1
			else:
				f1.close()
				cnt_file += 1
				cnt_img = 1
				f1 = open('/scratch365/sbanerj1/PaSC_2DAligned/pasc_control'+str(cnt_file)+'.txt','w')
				f1.write(file1+'\n')
				cnt_img += 1
				c += 1
print c
f1.close()
'''
