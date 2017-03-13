import os
import numpy as np
#import matplotlib.pyplot as plt
import sys

fold = '/scratch365/sbanerj1/PaSC_Old_Vito_Asym/Original'
f1 = open('/scratch365/sbanerj1/PaSC_Old_Vito_Asym/Original/Full_handheld_Feat_Orig.txt','w')

for i in range(1,11):
	strng = 'pasc_handheld_Feat' + str(i) + '.txt'
	f2 = open(strng,'r')
	print 'Opened ', strng
	for line in f2:
		f1.write(line)
	f2.close()
	
f1.close()