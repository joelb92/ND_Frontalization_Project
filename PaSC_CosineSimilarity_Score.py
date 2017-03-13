import numpy as np
import os
import sys
from scipy import spatial

f3 = open('C:/Users/Sandipan/Desktop/NDFaceNet_New/PaSC_Old_Vito/CMC_Control_Front.txt','a')

#f1 = open('/scratch365/sbanerj1/PaSC_2DAligned/PaSC_handheld_Feat/Full_handheld_Feat_Averaged.txt','r')
f1 = open('Full_handheld_Feat_Orig_Averaged.txt','r')
#f2 = open('/scratch365/sbanerj1/PaSC_2DAligned/PaSC_handheld_Feat/Full_handheld_CosSim1.txt','w')
f2 = open('Full_handheld_CosSim_sym.txt','w')

name_list = []
feat_list = []
name = 'Orig'

for line in f1:
	fl = line.split(',')
	name_list.append(fl[0])
	temp_list = []
	for i in range(1,len(fl)):
		temp_list.append(float(fl[i]))
	feat_list.append(temp_list)
	
def score_C(temp_s,cls,score_F,ht):
	flag = 0
	for k in range(0,ht):
		cls1 = temp_s[k].split('d')[0]
		if cls == cls1 and flag == 0:
			score_F[ht-1] += 1
			flag = 1
	return score_F

score1 = 0
score5 = 0
score10 = 0
score_F = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
vid = 1145.0
cnt = 0
for i in range(0,len(name_list)):
	cls = name_list[i].split('d')[0]
	cos_score = []
	f2.write(name_list[i])
	#f2.write(',')
	for j in range(0,len(name_list)):
		sim = 1 - spatial.distance.cosine(feat_list[i], feat_list[j])
		cos_score.append(sim)
		f2.write(',')
		f2.write(str(sim))
	f2.write('\n')
	s = sorted(range(len(cos_score)), key=lambda k: cos_score[k])
	s.reverse()
	flag = 0
	temp_s = []
	
	#score_F = score_C(temp_s,cls,score_F,ht)
	
	for k in range(1,16):
		#print 'k ', k
		#print name_list[i], name_list[s[k]]
		temp_s.append(name_list[s[k]])
		
	for k in range(1,16):
		score_F = score_C(temp_s,cls,score_F,k)
	cnt += 1
	print cnt

f3.write(name)
for i in range(0,len(score_F)):
	score_F[i] = float("{0:.4f}".format(score_F[i]/vid))
	f3.write(',')
	f3.write(str(score_F[i]*100))
f3.write('\n')
print 'Final score ', score_F
print score_F[0], score_F[4], score_F[9], score_F[14]

f1.close()
f2.close()

