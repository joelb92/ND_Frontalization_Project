# Python wrapper for dlib face landmark detector (Kazemi and Sullivan, CVPR'14)
import dlib
#from skimage import io
import numpy as np
#from Utils import HOME
import os
import sys
#import oc2 as io

def _shape_to_np(shape):
    xy = []
    for i in range(68):
        xy.append((shape.part(i).x, shape.part(i).y,))
    xy = np.asarray(xy, dtype='float32')
    return xy

def get_faces(img1d,w,h,detector):
    # predictor_path = 'shape_predictor_68_face_landmarks.dat' # http://sourceforge.net/projects/dclib/files/dlib/v18.10/shape_predictor_68_face_landmarks.dat.bz2
    # detector = dlib.get_frontal_face_detector()
    # predictor = dlib.shape_predictor(predictor_path)
    img = np.uint8(np.asarray(np.reshape(img1d, (int(h), int(w)))))
    # result = Image.fromarray(img.astype(np.uint8))

    # return type(img)
    # cv.imshow("test",img)
    # lmarks = []
    # bboxes = []
    # for i,line in enumerate(line_arr):
    #    print('%d/%d'%(i,len(line_arr)))
    dets = detector(img, 0)
    flatDets = []
    for k, d in enumerate(dets):
        flatDets.append(d.left())
        flatDets.append(d.top())
        flatDets.append(d.right()-d.left())
        flatDets.append(d.bottom()-d.top())
        # print("Detection {}: Left: {} Top: {} Right: {} Bottom: {}".format(
        #     k, d.left(), d.top(), d.right(), d.bottom()))
    xy_flat = np.int64(np.asarray(flatDets)).tolist()
    return xy_flat


def get_landmarks(img1d,w,h,detector,predictor):

    # predictor_path = 'shape_predictor_68_face_landmarks.dat' # http://sourceforge.net/projects/dclib/files/dlib/v18.10/shape_predictor_68_face_landmarks.dat.bz2
    # detector = dlib.get_frontal_face_detector()
    # predictor = dlib.shape_predictor(predictor_path)
    print('20')
    img = np.uint8(np.asarray(np.reshape(img1d,[int(h),int(w)])))
    # result = Image.fromarray(img.astype(np.uint8))

    # return type(img)
    # cv.imshow("test",img)
    # lmarks = []
    # bboxes = []
    #for i,line in enumerate(line_arr):
    #    print('%d/%d'%(i,len(line_arr)))
    dets = detector(img, 0)
#     print "found " + str(len(dets)) + " faces"
    if len(dets) == 0:
        rect = dlib.rectangle(0,0,img.shape[0], img.shape[1])
    else:
        rect = dets[0]

    	shape = predictor(img, rect)
    	xy = _shape_to_np(shape)
    	xy_flat = np.int64(np.reshape(xy,[1,68*2])).tolist()
    	np.savetxt('landmarks.txt', xy)
    	return xy_flat
    return np.int64(np.asarray([0,0])).tolist()
    # # lmarks.append(xy)
    # # bboxes.append(rect)
    #
    # lmarks = np.vstack(lmarks)
    # bboxes = np.asarray(bboxes)
    # print lmarks
    # return lmarks

if __name__ == "__main__":
    get_landmarks(sys.argv[1])
