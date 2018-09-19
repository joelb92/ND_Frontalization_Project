% DEMO script for frontalizing faces, based on the method described in [1] and [2]
% This demo produces a frontalized face for the image test.jpg
%
% If you find this code useful and use it in your own work, please add
% reference to [2] and, if appropriate, any of the other papers mentioned below.
%
% Dependencies:
%   The demo uses the following dependencies. You MUST have these installed and
%   available on the MATLAB path:
%   1. calib function available from [3,4], required for estimating the query
%       camera projection matrix, C_Q. Calib functions are assumed to be
%       under folder <frontalization>/calib
%   2. Facial feature detection functions. The demo provides examples
%       showing frontalization results obtained using the SDM method [5]
%       (default used in the paper), the facial feature detector of Zhu and
%       Ramanan [6], or the dlib detector (Kazemi and Sullivan) [7]. Please see
%       the script facial_feature_detection.m on how to use these (as well as
%       edit paths to the detector used in practice, in case these differ from
%       the ones in the script). See the function makeNew3DModel.m in case
%       a different facial feature detector is used.
%   3. OpenCV required by calib for calibration routines and some of the
%       detectors for cascase classifiers (e.g., SDM)
%	4. Dlib installed and added to your system PATH (can be installed using pip install dlib)
%
%  References:
%	[1] Banerjee, S., Brogan, J., Krizaj, J., Bharati, A., Webster, B.R., Struc, V., Flynn, P.J. and Scheirer, W.J.
%	2018, March. "To Frontalize or Not to Frontalize: Do We Really Need Elaborate Pre-processing 
%	to Improve Face Recognition?" IEEE Winter Conference on Applications of Computer Vision (WACV),
%	Page 20-29, Mar. 2018.
%
%   [2] Tal Hassner, Shai Harel, Eran Paz, Roee Enbar, "Effective Face
%   Frontalization in Unconstrained Images," forthcoming.
%   See project page for more details:
%   http://www.openu.ac.il/home/hassner/projects/frontalize
%
%   [3] T. Hassner, L. Assif, and L. Wolf, "When Standard RANSAC is Not Enough: Cross-Media
%   Visual Matching with Hypothesis Relevancy," Machine Vision and Applications (MVAP),
%   Volume 25, Issue 4, Page 971-983, 2014
%   Available: http://www.openu.ac.il/home/hassner/projects/poses/
%
%   [4] T. Hassner, "Viewing Real-World Faces in 3D," International Conference on Computer Vision (ICCV),
%   Sydney, Austraila, Dec. 2013
%   Available: http://www.openu.ac.il/home/hassner/projects/poses/
%
%   [5] X. Xiong and F. De la Torre, "Supervised Descent Method and its Application to Face
%   Alignment," IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2013
%   Available: http://www.humansensing.cs.cmu.edu/intraface
%
%   [6] X. Zhu, D. Ramanan. "Face detection, pose estimation and landmark localization in the wild,"
%   Computer Vision and Pattern Recognition (CVPR) Providence, Rhode Island, June 2012.
%   Available: http://www.ics.uci.edu/~xzhu/face/
%
%   [7] V. Kazemi, J. Sullivan. "One Millisecond Face Alignment with an
%   Ensemble of Regression Trees," Computer Vision and Pattern Recognition
%   (CVPR), Columbus, Ohio, June, 2014
%   Available through the dlib library: http://blog.dlib.net/2014/08/real-time-face-pose-estimation.html

%   The SOFTWARE ("frontalization" and all included files) is provided "as is", without any
%   guarantee made as to its suitability or fitness for any particular use.
%   It may contain bugs, so use of this tool is at your own risk.
%   We take no responsibility for any damage that may unintentionally be caused
%   through its use.
%
%   ver 1.2, 18-May-2015
%
% clear;
function Frontalize_Single(imagePath,output_dir,detector,frontalization,faceSize)
% close all;
% clear all;

addpath calib;
addpath parfor_progress;
addpath kakearney-subdir/subdir;
load eyemask eyemask
predictor_file = fullfile(pwd,'shape_predictor_68_face_landmarks.dat')
% dlibpredictor = py.dlib.shape_predictorf(predictor_file);
% dlibdetector = py.dlib.get_frontal_face_detector()
% mod = py.importlib.import_module('dlib_detect_script_optimized');
% mod = py.reload(mod);
% clusterName = parallel.importProfile(cluster_profile);


% Load query image
eyemask = load('eyemask');
eyemask = eyemask.eyemask;
% load eyemask eyemask % mask to exclude eyes from symmetry
% load DataAlign2LFWa REFSZ REFTFORM % similarity transf. from rendered view to LFW-a coordinates
Model3D = [];
if strcmp(detector,'ZhuRamanan')
    Model3D = load('model3DZhuRamanan.mat');
    Model3D = Model3D.Model3D;
elseif strcmp(detector,'dlib')
    Model3D = load('model3Ddlib.mat');
    Model3D = Model3D.model_dlib;
elseif strcmp(detector,'CMR')
    Model3D = load('model3DCMRfrontal.mat');
    Model3D = Model3D.Model3D;
end


others = load('DataAlign2LFWa');
REFSZ = others.REFSZ;
REFTFORM = others.REFTFORM;


% Model3D = load(strcat('model3D',detector,'.mat'));
% Model3D = Model3D.Model3D;
xcordsMask = double(Model3D.ref_XY(:,1));
ycordsMask = double(Model3D.ref_XY(:,2));
hull = convhull(xcordsMask,ycordsMask);
xcords = xcordsMask(hull);
ycords = ycordsMask(hull);
mask = poly2mask(xcords,ycords,REFSZ(1),REFSZ(2));
mask = imtransform(mask,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);
[maski,maskj] = find(mask);
maskBB = [min(maskj),min(maski),max(maskj)-min(maskj),max(maski)-min(maski)];
mask = uint8(reshape([mask mask mask],[size(mask) 3]));


tic
failedImgs = 0;
successImgs = 0;
numFinished = 0;

% load model3DZhuRamanan Model3D

if 1    
    numFinished = numFinished+1;

    
    
    
    filename = imagePath
    [full_dir,fname,ext] = fileparts(filename);
            base_dir = full_dir;
            remaining_dir_fail = full_dir(length(base_dir)+1:length(full_dir));
            save_dir_failed = fullfile(output_dir,'failed',[fname ext]);
            [faildir, ~, ~] = fileparts(save_dir_failed);
            mkdir(faildir);
            
    try
        I_Q = imread(filename);
    catch
        display('Image was in CMYK, skipping');
        failedImgs = failedImgs+1;
        try
            copyfile(filename,fullfile(save_dir_failed,strcat(fname,ext)),'f')
        catch
        end
        continue;
    end
    imSize = size(I_Q);
    if length(imSize) == 3
        if imSize(3) ~= 3
            I_Q=reshape([I_Q I_Q I_Q],[size(I_Q) 3]);
        end
    else
        I_Q=reshape([I_Q I_Q I_Q],[size(I_Q) 3]);
    end
    
    
    if isempty(I_Q)
        display('could not load image')
        failedImgs = failedImgs+1;
        try
            copyfile(filename,fullfile(save_dir_failed,strcat(fname,ext)),'f')
        catch
        end
        continue
    end
    
    
    
    switch frontalization,
        case 'Hassner',
            fidu_XY = [];
            frontal_raw = [];
            [full_dir,fname,ext] = fileparts(filename);
            base_dir = full_dir;
            remaining_dir = full_dir(length(base_dir)+1:length(full_dir));
            display(remaining_dir);
            save_dir_sym = fullfile(output_dir,'sym',[fname ext]);
            save_dir_asym = fullfile(output_dir,'asym',[fname ext]);
            save_dir_failed = fullfile(output_dir,'failed',[fname ext]);
            try
            [fidu_XY,frontal_raw,frontal_sym,hardsym_images]=demo(I_Q,detector,filename,[],faceSize,Model3D,'','','',eyemask);
            catch
                disp('failed')
                [full_dir,fname,ext] = fileparts(save_dir_failed);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
            try
                imwrite(I_Q,save_dir_failed)
            catch
            end
            end
            if ~isempty(fidu_XY) && ~isempty(frontal_raw)
                display ('Found face in image!')
            %         figure; imshow(I_Q); hold on; plot(fidu_XY(:,1),fidu_XY(:,2),'.');
            %         hold off; title('Query photo with detections overlaid');
            %         figure; imshow(frontal_raw); title('Frontalized no symmetry');
            %         figure; imshow(frontal_sym); title('Frontalized with soft symmetry');
            landmarks = fidu_XY;
            
            
           
            [full_dir,fname,ext] = fileparts(save_dir_sym);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
            [full_dir,fname,ext] = fileparts(save_dir_asym);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
            
            imwrite(frontal_raw,save_dir_asym);
            imwrite(frontal_sym,save_dir_sym);
            display(filename);
            display(save_dir_asym);
            else
                display('Could not detect face in image');
                failedImgs = failedImgs+1
                [full_dir,fname,ext] = fileparts(save_dir_failed);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
                 try
                     display(save_dir_failed)
                imwrite(I_Q,save_dir_failed);
                catch
                end
            end
            %         display(fidu_XY)
        case 'Vito',
%             addpath('Vito') ;
                 [full_dir,fname,ext] = fileparts(filename);
                base_dir = full_dir;
                remaining_dir = full_dir(length(base_dir)+1:length(full_dir));
                save_dir_asym = fullfile(output_dir,'asym',[fname ext]);
                save_dir_sym1 = fullfile(output_dir,'sym1',[fname ext]);
                save_dir_sym2 = fullfile(output_dir,'sym2',[fname ext]);
                save_dir_sym3 = fullfile(output_dir,'sym3',[fname ext]);
                save_dir_failed = fullfile(output_dir,'failed',[fname ext]);
                frontal_raw = [];
%                 try
            [warped_surface, fidu_XY,landmarked_img,frontal_raw,hardsym_images]=demoVito(I_Q,detector,filename,faceSize,Model3D,'','','');
%                 catch
                     disp('failed')
                     [full_dir,fname,ext] = fileparts(save_dir_failed);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
%             try
%                 imwrite(I_Q,save_dir_failed)
%             catch
%             end
%                 end
%             rawsize = length(frontal_raw)
%             display('made it here')
            if ~isempty(frontal_raw)
                display('frontalization succeeded!')
                X=landmarked_img.img;
                landmarks=landmarked_img.landmarks;
                bbox=landmarked_img.bbox;


                [full_dir,fname,ext] = fileparts(save_dir_sym1);
                
                if exist(full_dir, 'dir') == 0
                    mkdir(full_dir);
                end
                [full_dir,fname,ext] = fileparts(save_dir_sym2);
                
                if exist(full_dir, 'dir') == 0
                    mkdir(full_dir);
                end
                [full_dir,fname,ext] = fileparts(save_dir_sym3);
                
                if exist(full_dir, 'dir') == 0
                    mkdir(full_dir);
                end
                [full_dir,fname,ext] = fileparts(save_dir_asym);
                if exist(full_dir, 'dir') == 0
                    mkdir(full_dir);
                end
                            [full_dir,fname,ext] = fileparts(save_dir_failed);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
                imwrite(frontal_raw,save_dir_asym);
                imwrite(hardsym_images{1},save_dir_sym1);
                imwrite(hardsym_images{2},save_dir_sym2);
                imwrite(hardsym_images{3},save_dir_sym3);
                %         figure
                %         imshow(X,[])
                %         hold on
                % %         rectangle('Position',bbox,'EdgeColor','r','LineWidth',2);
                %         plot(landmarks(:,1)',landmarks(:,2)','go',...
                %             'MarkerSize',2,...
                %             'MarkerEdgeColor','g',...
                %             'MarkerFaceColor','g')
                %         display(fidu_XY)
                
%                 surf(warped_surface.X, warped_surface.Y, warped_surface.Z, warped_surface.C, 'edgecolor', 'none', 'FaceColor','texturemap')
            else
                display('Could not detect face in image: Saving');
                failedImgs = failedImgs+1
                 [full_dir,fname,ext] = fileparts(save_dir_failed);
            if exist(full_dir, 'dir') == 0
                display(strcat('Making dir: ', full_dir));
                mkdir(full_dir);
            end
                 try
                imwrite(I_Q,save_dir_failed);
                 catch
                 end
            end
        otherwise,
            disp('A combination of this Landmarker and Frontalization technique does not exist !!') ;
            
    end
    
    [full_dir,fname,ext] = fileparts(filename);
    base_dir = full_dir;
    remaining_dir = full_dir(length(base_dir)+1:length(full_dir));
    save_dir_sym = fullfile(output_dir,'sym',[fname ext]);
    save_dir_asym = fullfile(output_dir,'asym',[fname ext]);
    save_dir_failed = fullfile(output_dir,'failed',[fname ext]);
    [full_dir,fname,ext] = fileparts(save_dir_sym);
    
    if exist(full_dir, 'dir') == 0
        mkdir(full_dir);
    end
    [full_dir,fname,ext] = fileparts(save_dir_asym);
    if exist(full_dir, 'dir') == 0
        mkdir(full_dir);
    end
    
    
    
    
end
totalTime = toc;
disp('Finished all files!')
disp('Total time:')
disp(totalTime)
disp('Total images frontalized:')
disp(successImgs)
disp('Total images failed to frontalize:')
disp(failedImgs)

