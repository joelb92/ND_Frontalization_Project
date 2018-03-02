% Script used to demonstrate the use of different facial feature detectors
% for face frontalizations.
% PLEASE MODIFY THIS SCRIPT AS REQUIRED BY THE FACIAL FEATURE DETECTOR
% INSTALLED ON YOUR SYSTEM. This is part of the distribution for
% face image frontalization ("frontalization" software), described in [1].
%
% If you find this code useful and use it in your own work, please add
% reference to [1]. Please also respect any distribution notices of the
% facial feature detectors used [2,3,4]
%
% Please see project page for more details:
%   http://www.openu.ac.il/home/hassner/projects/frontalize
%
% Please see demo.m for example usage.
%
%  References:
%   [1] Tal Hassner, Shai Harel, Eran Paz, Roee Enbar, "Effective Face
%   Frontalization in Unconstrained Images," forthcoming.
%   See project page for more details:
%   http://www.openu.ac.il/home/hassner/projects/frontalize
%
%   [2] X. Xiong and F. De la Torre, "Supervised Descent Method and its Application to Face
%   Alignment," IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2013
%   Available: http://www.humansensing.cs.cmu.edu/intraface
%
%   [3] X. Zhu, D. Ramanan. "Face detection, pose estimation and landmark localization in the wild,"
%   Computer Vision and Pattern Recognition (CVPR) Providence, Rhode Island, June 2012.
%   Available: http://www.ics.uci.edu/~xzhu/face/
%
%   [4] V. Kazemi, J. Sullivan. "One Millisecond Face Alignment with an
%   Ensemble of Regression Trees," Computer Vision and Pattern Recognition
%   (CVPR), Columbus, Ohio, June, 2014
%   Available through the dlib library: http://blog.dlib.net/2014/08/real-time-face-pose-estimation.html
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
function [Model3D,fidu_XY,bbox] = facial_feature_detection(detector,I_Q,filename,landmarks,Model3D,dlibdetector,dlibpredictor,dlibmodule,detMod)
fidu_XY = [];
mod = dlibmodule;
bbox = [];
switch detector
    case 'SDM'
        % SDM detector (Intraface) [2]
        if isempty(Model3D)
            load model3DSDM Model3D % reference 3D points corresponding to SDM detections
        end
        % detect facial features on query
        addpath(genpath('SDM/FacialFeatureDetection&Tracking_v1.4'))
        tmpdir = pwd;
        cd SDM/FacialFeatureDetection&Tracking_v1.4
        [Models,option] = xx_initialize;
        cd(tmpdir);
        faces = Models.DM{1}.fd_h.detect(I_Q,'MinNeighbors',option.min_neighbors,...
            'ScaleFactor',1.2,'MinSize',[50 50]);
        if ~isempty(faces)
            output = xx_track_detect(Models,I_Q,faces{1},option);
            fidu_XY = double(output.pred);
        end
        
        
        
    case 'ZhuRamanan'
        % Zhu and Ramanan detector [3]
        %if isempty(Model3D)
            load model3DZhuRamanan Model3D % reference 3D points corresponding to Zhu & Ramanan detections
        %end
%         Model3D = []
        if ~isempty(landmarks)
            fidu_XY = landmarks
            return
        end
        % detect facial features on query
        addpath(genpath('ZhuRamanan'))
        model =[]
        if isempty(model)
        load('ZhuRamanan/face_p146_small.mat','model');
        end
        model.interval = 5;
        model.thresh = min(-0.65, model.thresh);
        if length(model.components)==13
            posemap = 90:-15:-90;
        elseif length(model.components)==18
            posemap = [90:-15:15 0 0 0 0 0 0 -15:-15:-90];
        else
            error('Can not recognize this model');
        end
        
        
        I_Q_bs = detect(I_Q, model, model.thresh);
        if isempty(I_Q_bs)
            return
        end
        
        I_Q_bs = clipboxes(I_Q, I_Q_bs);
        I_Q_bs = nms_face(I_Q_bs,0.3);
        
        if (isempty(I_Q_bs))
            return;
        end
        x1 = I_Q_bs(1).xy(:,1);
        y1 = I_Q_bs(1).xy(:,2);
        x2 = I_Q_bs(1).xy(:,3);
        y2 = I_Q_bs(1).xy(:,4);
        fidu_XY = [(x1+x2)/2,(y1+y2)/2];
        fidu_XY = uint32(fidu_XY);
        minx = min(fidu_XY(1,:));
        miny = min(fidu_XY(1,:));
        bw = max(fidu_XY(1,:))-minx;
        bh = max(fidu_XY(2,:))-miny;
        bbox = [minx,miny,bw,bh];
        %         zhurama_hard_symmetry(I_Q,fidu_XY,Model3D);
        
    case 'dlib'
        % Kazemi and Sullivan detector implemented by the dlib library [4]
        if isempty(Model3D)
            load model3Ddlib1 Model3D % reference 3D points corresponding to dlib detections
        end
        
        %addpath('/afs/nd.edu/user25/jbrogan4/Public/Software/dlib/')
        
        % this detector runs in Python; here we only load them for the
        % current face image. See dlib_detect_script.py for an example
        % usage
        
        %% commands to call the python script for dlib landmark detection from matlab
        %         %% does not work as Utils and scipy not compiled on server
        predictor_file = fullfile(pwd,'shape_predictor_68_face_landmarks.dat');
        predictor = py.dlib.shape_predictor(predictor_file);
        detector = py.dlib.get_frontal_face_detector();
        mod = py.importlib.import_module('dlib_detect_script_optimized');
        mod = py.reload(mod);
        %          detector = dlibdetector;
        %          predictor = dlibpredictor;
        gimg = uint8(rgb2gray(I_Q));
        
        imSize = size(gimg);
        gimg = reshape(gimg',[1,imSize(1)*imSize(2)]);
%         img = np.uint8(np.asarray(np.reshape(img1d,[h,w])))
        points = mod.get_landmarks(gimg,int32(imSize(2)),int32(imSize(1)),detector,predictor);
        rawPoints = [];
        pointlength = length(points{1});
        if length(points{1}) > 2
        for i = 1:length(points{1})
            rawPoints = [rawPoints points{1}{i}];
        end
        fidu_XY = reshape(rawPoints,[2,length(points{1})/2])';
       else
            fidu_XY = [];
        end
        
        %         commandStr = ['python dlib_detect_script.py ',filename];
        %         [status, commandOut] = system(commandStr);
        %         status;
        %         if status==1
        %             fprintf(commandOut);
        %         end
        % %         fprintf(commandOut);
        % %       fidu_XY = load('dlib_xy.mat'); % load detections performed by Python script on current image
        %
        % %          [points,~]=dlibDetectFaceLandmarks(rgb2gray(I_Q),[],1);
        % %          fidu_XY = points;
        %         fidu_XY=[];
        %         fid = fopen('landmarks.txt','r');
        %         tline = fgetl(fid);
        %         while ischar(tline)
        %             strarr=strsplit(tline);
        %             if length(strarr) >= 2
        %             arr=[str2num(strarr{1,1}) str2num(strarr{1,2})];
        %             fidu_XY=[fidu_XY;arr];
        %             end
        %             tline = fgetl(fid);
        %         end
        %
        %         fclose(fid);
        if length(fidu_XY) == 68
            fidu_XY = reshape(fidu_XY,68,2);
            fidu_XY = double(fidu_XY);
            minx = min(fidu_XY(1,:));
            miny = min(fidu_XY(1,:));
            bw = max(fidu_XY(1,:))-minx;
            bh = max(fidu_XY(2,:))-miny;
            bbox = [minx,miny,bw,bh];
            %         zhurama_hard_symmetry(I_Q,fidu_XY,Model3Dw);
        end
    case 'CMR'
        % do the detection
        if isempty(Model3D)
            load model3DCMRfrontal Model3D
        end
        addpath calib
        addpath Vito
        mpath = 'model/';
        model = generateCMR_model1(mpath);
        FDetect = vision.CascadeObjectDetector;
        BB = step(FDetect,I_Q);
%         bbox = detect_faces_fp(I_Q,model.VJdetector,1);
        faces = BB
        if size(BB,1) > 0
        bbox=BB(1,:);
        
        % find the landmarks
        landmarks = CMRfind_facial_landmarks(I_Q,bbox,model.Ldetector,1);
        fidu_XY=reshape(landmarks,[68,2]);
        % rescale image and landmarks
        bbox_width = bbox(3);
        bbox_norm = 300;
        scale_fac = bbox_norm/bbox_width;
        bbox=bbox*scale_fac;
        %         fidu_XY=landmarks*scale_fac;
        %         X=imresize(X,scale_fac,'bilinear');
        else
            display('Could not find any faces!')
            bbox = []
        end
    otherwise
        error(1,'To use a new, unsupported facial feature detector please see MakeNew3DModel.m\n');
        
end