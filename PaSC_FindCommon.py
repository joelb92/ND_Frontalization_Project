import os
import sys
from shutil import copytree

f1 = open('Full_control_Feat_Pre_Averaged.txt','r')
f2 = open('pasc_intersection.txt','r')
f3 = open('pasc_control_feat_Pre.txt','w')

name1 = []
name2 = []
line1 = []

for line in f1:
	fl = line.split(',')
	name1.append(fl[0])
	line1.append(line)
for line in f2:
	fl = line.split('\n')
	name2.append(fl[0])
cnt = 0
for i in range(0,len(name1)):
	if name1[i] in name2:
		f3.write(line1[i])
		cnt += 1
print cnt