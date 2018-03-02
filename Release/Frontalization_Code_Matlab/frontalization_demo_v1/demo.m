% DEMO script for frontalizing faces, based on the method described in [1]
% This demo produces a frontalized face for the image test.jpg
%
% If you find this code useful and use it in your own work, please add
% reference to [1] and, if appropriate, any of the other papers mentioned below.
%
% Dependencies: 
%   The demo uses the following dependencies. You MUST have these installed and
%   available on the MATLAB path:
%   1. calib function available from [2,3], required for estimating the query
%       camera projection matrix, C_Q. Calib functions are assumed to be
%       under folder <frontalization>/calib
%   2. Facial feature detection functions. The demo provides examples
%       showing frontalization results obtained using the SDM method [4] 
%       (default used in the paper), the facial feature detector of Zhu and 
%       Ramanan [5], or the dlib detector (Kazemi and Sullivan) [6]. Please see
%       the script facial_feature_detection.m on how to use these (as well as 
%       edit paths to the detector used in practice, in case these differ from 
%       the ones in the script). See the function makeNew3DModel.m in case
%       a different facial feature detector is used.
%   3. OpenCV required by calib for calibration routines and some of the
%       detectors for cascase classifiers (e.g., SDM)
%
%  References:
%   [1] Tal Hassner, Shai Harel, Eran Paz, Roee Enbar, "Effective Face
%   Frontalization in Unconstrained Images," forthcoming. 
%   See project page for more details: 
%   http://www.openu.ac.il/home/hassner/projects/frontalize
%
%   [2] T. Hassner, L. Assif, and L. Wolf, "When Standard RANSAC is Not Enough: Cross-Media 
%   Visual Matching with Hypothesis Relevancy," Machine Vision and Applications (MVAP), 
%   Volume 25, Issue 4, Page 971-983, 2014 
%   Available: http://www.openu.ac.il/home/hassner/projects/poses/
%
%   [3] T. Hassner, "Viewing Real-World Faces in 3D," International Conference on Computer Vision (ICCV), 
%   Sydney, Austraila, Dec. 2013
%   Available: http://www.openu.ac.il/home/hassner/projects/poses/
%
%   [4] X. Xiong and F. De la Torre, "Supervised Descent Method and its Application to Face
%   Alignment," IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2013
%   Available: http://www.humansensing.cs.cmu.edu/intraface
%
%   [5] X. Zhu, D. Ramanan. "Face detection, pose estimation and landmark localization in the wild," 
%   Computer Vision and Pattern Recognition (CVPR) Providence, Rhode Island, June 2012. 
%   Available: http://www.ics.uci.edu/~xzhu/face/
%
%   [6] V. Kazemi, J. Sullivan. "One Millisecond Face Alignment with an
%   Ensemble of Regression Trees," Computer Vision and Pattern Recognition
%   (CVPR), Columbus, Ohio, June, 2014
%   Available through the dlib library: http://blog.dlib.net/2014/08/real-time-face-pose-estimation.html
%
%   Copyright 2014, Tal Hassner
%   http://www.openu.ac.il/home/hassner/projects/frontalize
%
%
%   The SOFTWARE ("frontalization" and all included files) is provided "as is", without any
%   guarantee made as to its suitability or fitness for any particular use.
%   It may contain bugs, so use of this tool is at your own risk.
%   We take no responsibility for any damage that may unintentionally be caused
%   through its use.
%
%   ver 1.2, 18-May-2015
%
% clear;
% close all;

function [fidu_XY,frontal_raw,frontal_sym,hardsym_images]=demo(I_Q,landmarker,file_path,landmark_detections,faceSize,Model3D,dlibdetector,dlibpredictor,dlibmodule,eyemask)
addpath calib
fidu_XY = [];
frontal_raw = [];
frontal_sym = [];
hardsym_images = {};
% Load query image
%I_Q =imresize(imread('/Users/abharati/Desktop/1st Sem/dlib-18.17/examples/build1/JPEG/02463d3207.jpg'),[250 250]);
% load some data
% load eyemask eyemask % mask to exclude eyes from symmetry
load DataAlign2LFWa REFSZ REFTFORM % similarity transf. from rendered view to LFW-a coordinates

% Detect facial features with prefered facial feature detector 
detector = landmarker; % alternatively 'ZhuRamanan', 'dlib'
% Note that the results in the paper were produced using SDM. We have found
% other detectors to produce inferior frontalization results. 
fidu_XY = [];
[Model3D,fidu_XY,bbox] = facial_feature_detection(detector,I_Q,file_path,landmark_detections,[],'','',dlibmodule);
xcordsMask = double(Model3D.ref_XY(:,1));
ycordsMask = double(Model3D.ref_XY(:,2));
hull = convhull(xcordsMask,ycordsMask);
xcords = xcordsMask(hull);
ycords = ycordsMask(hull);
mask = poly2mask(xcords,ycords,REFSZ(1),REFSZ(2));
mask = imtransform(mask,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);
[maski,maskj] = find(mask);
maskBB = [min(maskj),min(maski),max(maskj)-min(maskj),max(maski)-min(maski)];
if length(size(I_Q)) == 3
mask = uint8(reshape([mask mask mask],[size(mask) 3]));
end

if isempty(fidu_XY)
    display('Failed to detect facial features / find face in image.');
    lsize = size(fidu_XY)
    
else
%     imshow(I_Q)
%     hold on
%     plot(fidu_XY(:,1),fidu_XY(:,2),'.');
%     w = waitforbuttonpress;
% [I_Q roi]= maskFaceByPoints(I_Q,fidu_XY);
% Estimate projection matrix C_Q
failed = 0;
C_Q = [];
try
 [C_Q, ~,~,~] = estimateCamera(Model3D, fidu_XY);   
catch
    failed = 1
end
   if failed
       display('estimateCamera failed')
frontal_raw = [];
frontal_sym = [];

hardsym_images = {[],[],[],[]};
   else
       display('estimateCamera succeeded!')
        
%C_Q=[520,30,30,85810; 60,220,-470,118570; 0,0,0,770];

% Render frontal view
[frontal_sym, frontal_raw] = Frontalize(C_Q, I_Q, Model3D.refU, eyemask);


% Apply similarity transform to LFW-a coordinate system, for compatability
% with existing methods and results
frontal_sym = imtransform(frontal_sym,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);
frontal_raw = imtransform(frontal_raw,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);
frontal_sym = mask.*frontal_sym;
frontal_sym = imcrop(frontal_sym,maskBB);
frontal_raw = mask.*frontal_raw;
frontal_raw = imcrop(frontal_raw,maskBB);

flipped_raw=fliplr(frontal_raw);
flipped_sym=fliplr(frontal_sym);

hardsym_images = cell(4);
maped1 = [flipped_raw(:,1:size(frontal_raw,2)/2,:), frontal_raw(:,size(frontal_raw,2)/2+1:end,:)];
maped1 = uint8(maped1);
maped2 = [frontal_raw(:,1:size(frontal_raw,2)/2,:), flipped_raw(:,size(frontal_raw,2)/2+1:end,:)];
maped2 = uint8(maped2);
maped3 = [flipped_sym(:,1:size(frontal_sym,2)/2,:), frontal_sym(:,size(frontal_sym,2)/2+1:end,:)];
maped3 = uint8(maped3);
maped4 = [frontal_sym(:,1:size(frontal_sym,2)/2,:), flipped_sym(:,size(frontal_sym,2)/2+1:end,:)];
maped4 = uint8(maped4);

maped1 = imresize(maped1,faceSize);
maped2 = imresize(maped2,faceSize);
maped3 = imresize(maped3,faceSize);
maped4 = imresize(maped4,faceSize);

frontal_raw = imresize(frontal_raw,faceSize);
frontal_sym = imresize(frontal_sym,faceSize);

hardsym_images = {maped1,maped2,maped3,maped4};
   end
end
% figure
% subplot(2,4,1)
% imshow(uint8(I_Q),[])
% title('Original')
% 
% subplot(2,4,2)
% imshow(uint8(frontal_raw),[])
% title('Without Symmetry')
% 
% subplot(2,4,3)
% imshow(uint8(frontal_sym),[])
% title('With Soft Symmetry')
% 
% subplot(2,4,4)
% maped1 = [flipped_raw(:,1:size(frontal_raw,2)/2,:), frontal_raw(:,size(frontal_raw,2)/2:end,:)];
% imshow(uint8(maped1),[])
% title('With Hard Symmetry 1 (on image 2)')
% 
% subplot(2,4,5)
% maped2 = [frontal_raw(:,1:size(frontal_raw,2)/2,:), flipped_raw(:,size(frontal_raw,2)/2:end,:)];
% imshow(uint8(maped2),[])
% title('With Hard Symmetry 2 (on image 2)')
% 
% subplot(2,4,6)
% maped3 = [flipped_sym(:,1:size(frontal_raw,2)/2,:), frontal_sym(:,size(frontal_raw,2)/2:end,:)];
% imshow(uint8(maped3),[])
% title('With Hard Symmetry 1 (on image 3)')
% 
% subplot(2,4,7)
% maped4 = [frontal_sym(:,1:size(frontal_raw,2)/2,:), flipped_sym(:,size(frontal_raw,2)/2:end,:)];
% imshow(uint8(maped4),[])
% title('With Hard Symmetry 2 (on image 3)')
% Display results
% figure; imshow(I_Q); title('Query photo');
% figure; imshow(I_Q); hold on; plot(fidu_XY(:,1),fidu_XY(:,2),'.'); hold off; title('Query photo with detections overlaid');
% figure; imshow(frontal_raw); title('Frontalilzed no symmetry');
% figure; imshow(frontal_sym); title('Frontalilzed with soft symmetry');


end
