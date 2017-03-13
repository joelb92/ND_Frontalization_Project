import numpy as np
import os
import cv2
import sys
import time
import copy
import random
import math
from cv2 import cv
from scipy.misc import imread
#from matplotlib import pyplot as ppl

dlibf = '/home/sandipan/CW_CMR_Hass/Sym/sym'

traintxt = '/home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Train.txt'
val1txt = '/home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Val1.txt'
#val2txt = '/home/sandipan/SREFI_Full_Set/SREFI_Full_NoMask_Val2_1.txt'

trainf = '/home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Train'
val1f = '/home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Val1'
print 'blah New'
#val2f = '/home/sandipan/SREFI_Full_Set/SREFI_Full_NoMask_Val2_'

if not os.path.exists(trainf):
    os.makedirs(trainf)
if not os.path.exists(val1f):
    os.makedirs(val1f)
#if not os.path.exists(val2f):
#    os.makedirs(val2f)

f1 = open(traintxt,'w')
f2 = open(val1txt,'w')
#f3 = open(val2txt,'w')

label = 0
for root,dirs,files in os.walk(dlibf):
    for dirc in dirs:
    	#print 'dir ', dirc
        length = len(os.walk(dlibf+'/'+dirc).next()[2])
        #print 'length ', length
        arr = random.sample(range(1, (length+1)), length)
        train = int(math.ceil(0.9*length))
        #print train
        val = length
        #test = length
        arr1 = range(0, train)
        arr2 = range(train,length)
        #arr2 = arr[train:val]
        #print arr1, arr2
        #arr3 = arr[val:test]
        count = 0
        files1 = os.listdir(dlibf+'/'+dirc)
        for file1 in files1:
            if file1.endswith('.jpg') and cv2.imread(dlibf+ '/' + dirc + '/' + file1) is not None:
                strng = dlibf+ '/' + dirc + '/' + file1
                img = cv2.imread(strng)
                r_img = cv2.resize(img, (256, 256))
                if count in arr1:
                            strng0 = dirc + '_' + file1
                            strng1 = trainf + '/' + strng0
                            cv2.imwrite(strng1,r_img)
                            f1.write(strng0)
                            f1.write(' ')
                            f1.write(str(label))
                            f1.write('\n')
                elif count in arr2:
                            strng0 = dirc + '_' + file1
                            strng1 = val1f + '/' + strng0
                            cv2.imwrite(strng1,r_img)
                            f2.write(strng0)
                            f2.write(' ')
                            f2.write(str(label))
                            f2.write('\n')
#                 elif count in arr3:
#                             strng0 = dirc + '_' + file1
#                             strng1 = val2f + '/' + strng0
#                             cv2.imwrite(strng1,r_img)
#                             f3.write(strng0)
#                             f3.write(' ')
#                             f3.write(str(label))
#                             f3.write('\n')
                count += 1
            else:
                print 'Image not read'
        print 'count = ', count
        print 'label ', label
        label += 1
                
print 'label = ', label
f1.close()
f2.close()
#f3.close()                

        
        
