import os
import sys
from shutil import copytree

f1 = open('/afs/crc.nd.edu/user/s/sbanerj1/vol02/FaceNet_Project/Intersection/pasc_intersection.txt','r')
fold1 = sys.argv[1]
fold2 = sys.argv[2]
cnt = 0
for line in f1:
	name = line.split('\n')[0]
	src = fold1 + '/' + name
	if os.path.exists(src):
		dst = fold2 + '/' + name
		copytree(src,dst)
		cnt += 1

print cnt
