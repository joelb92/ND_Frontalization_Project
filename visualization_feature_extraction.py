import os
import numpy as np
#import matplotlib.pyplot as plt
import sys
import caffe

def main(argv):
	pasc_fold = argv[1]  # PaSC folder containing images
	#img_file = argv[2]   # text file containing image paths for feature extraction
	caffe_dep = argv[2]  # deploy.prototxt for caffe network
	caffe_mod = argv[3]  # .caffemodel file for trained network
	out_file = argv[4]   # full path of output file containing features
	ext_strng = argv[5]  # image extension ('.png' or '.jpg' or '.bmp' ...)
	
	filename1 = []
	filename2 = []
	
	fls = os.listdir(pasc_fold)
	for root,dirs,files in os.walk(pasc_fold):
		for dir in dirs:
			fls = os.listdir(pasc_fold+'/'+dir)
			for f in fls:
				if f.endswith(ext_strng):
					strng = pasc_fold+'/'+dir+'/'+f
					filename1.append(strng)
					filename2.append(f)
	# for fl1 in fls:
	# 	if fl1.endswith(ext_strng):
	# 		strng = pasc_fold + '/' + fl1
	# 		filename1.append(strng)
	# 		filename2.append(fl1)
			
	'''
	f1 = open(img_file,'r')
	filename1 = []
	filename2 = []
	for line in f1:
		fold = line.split('-')[0]
		img = line.split('\n')[0]
		strng = pasc_fold + '/' + fold + '/' + img
		filename1.append(strng)
		filename2.append(img)
	'''	
	caffe.set_mode_cpu()
	net = caffe.Net( caffe_dep, caffe_mod, caffe.TEST)
	transformer = caffe.io.Transformer({'data': net.blobs['data'].data.shape})
	transformer.set_transpose('data', (2,0,1))
	transformer.set_raw_scale('data', 255)  
	transformer.set_channel_swap('data', (2,1,0))
	img_dim = 224 #sub-crop dimensions used by model (224 for VGG, 227 for AlexNet/CaffeNet)
	net.blobs['data'].reshape(1,3,img_dim,img_dim)
	
	f2 = open(out_file,'w')
	
	for i in range(0,len(filename1)):
		if os.path.getsize(filename1[i]) > 0:
			net.blobs['data'].data[...] = transformer.preprocess('data', caffe.io.load_image(filename1[i]))
			out = net.forward()
			[(k, v.data.shape) for k, v in net.blobs.items()]
			[(k, v[0].data.shape) for k, v in net.params.items()]
			feat = net.blobs['fc7'].data[0]
			f2.write(filename2[i])
			f2.write(',')
			for tVal in range(0,len(feat)):
				f2.write(str(feat[tVal]))
				if tVal < (len(feat)-1):
					f2.write(',')
			f2.write('\n')
		print i
		
	f2.close()
	#f1.close()
	#print 'Feat Ext Done for ', 
main(sys.argv)