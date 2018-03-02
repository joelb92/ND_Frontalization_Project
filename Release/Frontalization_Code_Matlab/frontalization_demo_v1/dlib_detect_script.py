# Python wrapper for dlib face landmark detector (Kazemi and Sullivan, CVPR'14)
import dlib
from skimage import io
import glob
#from skimage import io
import numpy as np
#from Utils import HOME
from scipy.io.matlab import savemat
import os
import sys
#import oc2 as io

def _shape_to_np(shape):
    xy = []
    for i in range(68):
        xy.append((shape.part(i).x, shape.part(i).y,))
    xy = np.asarray(xy, dtype='float32')
    return xy


def get_landmarks(line):
    predictor_path = 'shape_predictor_68_face_landmarks.dat' # http://sourceforge.net/projects/dclib/files/dlib/v18.10/shape_predictor_68_face_landmarks.dat.bz2
    detector = dlib.get_frontal_face_detector()
    predictor = dlib.shape_predictor(predictor_path)

    lmarks = []
    bboxes = []
    #for i,line in enumerate(line_arr):
    #    print('%d/%d'%(i,len(line_arr)))
    img = io.imread(line)
    dets = detector(img, 0)
    if len(dets) == 0:
        rect = dlib.rectangle(0,0,img.shape[0], img.shape[1])
    else:
        rect = dets[0]

    shape = predictor(img, rect)
    xy = _shape_to_np(shape)
    np.savetxt('landmarks.txt', xy)
    
    lmarks.append(xy)
    bboxes.append(rect)

    lmarks = np.vstack(lmarks)
    bboxes = np.asarray(bboxes)
    print lmarks
    return lmarks

if __name__ == "__main__":
    get_landmarks(sys.argv[1])