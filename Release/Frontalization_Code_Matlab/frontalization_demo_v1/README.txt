Frontalization Demo V1
By Aparna Bharati and Joel Brogan

*****************************************
*Main function (Located in root folder)*
*****************************************
demo_aggregate(img_file,landmark_method,frontalization_method,verbosity)

landmark_method:
'ZhuRamanan' = Zhu and Ramanan landmarking ( mixtures of trees with a shared pool of parts by Xiangxin Zhu and Deva Ramanan [1]
'dlib' = Dlib landmarking (Ensemble of Regression Trees by Vahid Kazemi and Josephine Sullivan) [2]
'CMR' = CMR landmarking (Cascade Mixture of Regressors by Vito Struc) [NA]

frontalization_method:
'Hassner' = Effective Face Frontalization in Unconstrained Images by Tal Hassner et al. [3]
'Vito' = 2d to 3d camera estimation and Delaunay Triangulation [NA]

verbosity:
0 = view all possible frontalizations
1 = predict best frontalization based on Zhu Ramanan yaw estimation
	Predictions are as follows:
	0 to 15 degrees: Raw frontalized image (no symmetry used)
	15 to 45 degrees: Soft symmetry frontalized image
	45 to 90 degrees (if face is landmarkable): Hard symmetry frontalized image


*****************************************
*Images and Data						*
*****************************************
The images used to create the figures in the slides are held in the 'demoimages/' folder
Other usable images for the demo are held in the 'RFIIRimages/' folder

Images are named with the following rule: <ID><Pose>.jpg
ID = name
Pose = C, R, L for center, left, or right, respectively


*****************************************
*Important Demo Notes					*
*****************************************
-Currently, Zhu Ramanan + Vito frontalization is bugged, and outputting a garbled face.  It is suggested not to use this
-To reduce computation time, input images should first be cropped so most extraneous background noise is cut out of the image
-Using any 'Vito' based method will also produce a 3D face figure that can be rotated to see different views
-Using the prediction method will result in a slower calculation due to the added step of detecting yaw angle.  Currently, this step is performed using Zhu and Ramanan yaw estimation, which is slow.
-From our personal observations, when using Hassner's frontalization method, soft symmetry almost always seems to create a more plausible face


Sources:
[1]http://www.ics.uci.edu/~xzhu/face/
[2]http://sicv.activearchives.org/logbook/one-millisecond-face-alignment-with-an-ensemble-of-regression-trees/
[3]http://www.openu.ac.il/home/hassner/projects/frontalize/
