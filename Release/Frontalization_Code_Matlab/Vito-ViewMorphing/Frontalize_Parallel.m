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
function Frontalize_Parallel(image_dir,image_type,output_dir,detector,frontalization,numCores,faceSize,progressDir,jobName,numjobs,jobnum,shouldUpdateProgress)
% close all;
% clear all;

addpath parfor_progress;
addpath kakearney-subdir/subdir;
addpath('face_and_landmark_detector/')
% Load models
models = generateCMRmodel('models/');
predictor_file = fullfile(pwd,'shape_predictor_68_face_landmarks.dat')
% dlibpredictor = py.dlib.shape_predictor(predictor_file);
% dlibdetector = py.dlib.get_frontal_face_detector()
% mod = py.importlib.import_module('dlib_detect_script_optimized');
% mod = py.reload(mod);
% clusterName = parallel.importProfile(cluster_profile);
directoryPath = fullfile(image_dir);%'/Users/joel/Downloads/icbrw_Data/icbrw_ProbeImages/';
imagefiles = getLevel2FolderDifference(image_dir,fullfile(output_dir,'asym'),jobName,'.jpg',shouldUpdateProgress);
if ~strcmp(image_dir(end),'/')
    image_dir=strcat(image_dir,'/');
end
imagefiles = strcat(image_dir,imagefiles);
% for i=1:numel(imagefiles)
%     imagefiles{i}=fullfile(image_dir,imagefiles{i})
% end
jobnum = max(0,jobnum-1)
nfiles = numel(imagefiles);
stride = nfiles/numjobs;
fstart = ceil(1+stride*jobnum);
fend = ceil(min(1+fstart+stride,nfiles));
imagefiles=imagefiles(fstart:fend);
nfiles = numel(imagefiles);
display(strcat('image files total:',num2str(nfiles)))
display(strcat('running files from index ',num2str(fstart),' to ',num2str(fend)))
% reducedFile = strcat('fileNamesReduced_',jobName,'.mat');
% if ~exist(strcat('fileNames_',jobName,'.mat'), 'file')
%     disp('loading all filenames...')
%     folderLevel1 = dir(fullfile(directoryPath));
%     imagefiles = [];
%     
%     for n = 1:length(folderLevel1)
%         folder = folderLevel1(n);
%         if folder.isdir
%             imagefilestmp = dir(fullfile(directoryPath,folder.name,['*' image_type]));
%             for i = 1:length(imagefilestmp)
%                 imagefilestmp(i).name = fullfile(directoryPath,folder.name,imagefilestmp(i).name);
%             end
%             imagefiles = [imagefiles ; imagefilestmp];
%         end
%     end
%     
%     save(strcat('fileNames_',jobName),'imagefiles');
%     
% else
%     
%     if exist(reducedFile,'file') == 2
%         disp('loading filenames from reduced save file...')
%         load(reducedFile);
%     else
%         disp('loading filenames from save file...')
%         fileData = load(strcat('fileNames_',jobName));
%         fname = fieldnames(fileData(1));
%         imagefiles = getfield(fileData,fname{1})
%     end
% end
% 
% progressFileRegex = strcat('progress_',jobName,'*.txt');
% progressFiles = dir(progressFileRegex);
% alreadyProcessedIndexes = [];
% display('Modifying array to reflect progress from progress files...')
% for i = 1:length(progressFiles)
%     display(['found progress file ' progressFiles(i).name])
%     A = regexp( fileread(progressFiles(i).name), '\n', 'split');
%     if length(A) >= 2
%         start = str2num(A{2});
%         stop = str2num(A{1});
%         display(['Start: ' num2str(start) ' Stop: ' num2str(stop)])
%         alreadyProcessedIndexes = [alreadyProcessedIndexes [start:stop]];
%     end
% end
% if ~isempty(alreadyProcessedIndexes)
%     display('Saving reduced file list')
%     oldLength = length(imagefiles)
%     imagefiles = imagefiles(setdiff(1:length(imagefiles),alreadyProcessedIndexes));
%     newLength = length(imagefiles)
%     save(reducedFile,'imagefiles');
%     
% end
% delete(progressFileRegex);
% disp('fould all files...');
% % imagefiles = subdir(fullfile(directoryPath,['*' image_type]));
% disp(['done reading in ' int2str(size(imagefiles)) ' files!'])
% nfiles = length(imagefiles);    % Number of files found

% parpool('cvrl12',numCores)
% parpool('cvrl12',numCores)

% poolobj = gcp('nocreate');
predictor_file = fullfile(pwd,'shape_predictor_68_face_landmarks.dat');
dlibpredictor = py.dlib.shape_predictor(predictor_file);
dlibdetector = py.dlib.get_frontal_face_detector();
dlibmod = py.importlib.import_module('dlib_detect_script_optimized');
dlibmod = py.reload(dlibmod);
      

tic
failedImgs = 0;
successImgs = 0;
numFinished = 0;
parfor_progress(nfiles,progressDir);
 dmodel = load('ZhuRamanan/face_p146_small.mat');
 dmodel = dmodel.model;
% load('ZhuRamanan/face_p146_small.mat','model');
% parfor i=1:1
% end
% addAttachedFiles(gcp,'ZhuRamanan/face_p146_small.mat');

for i = 1:nfiles
    
%     percent = parfor_progress(-1,progressDir);
     numFinished = numFinished+1;
%     if 1
%         %keep track of progress
%         t = getCurrentTask();
%         progressFileName = strcat('progress_',jobName,'_',int2str(1),'.txt');
%         if exist(progressFileName,'file') == 2
%             A = regexp( fileread(progressFileName), '\n', 'split');
%             A{2} = sprintf('\n%i',i);
%             fid = fopen(progressFileName, 'w');
%             fprintf(fid, '%s', A{:});
%             fclose(fid);
%         else
%             fid = fopen(progressFileName, 'w');
%             fprintf(fid, '%i\n',i);
%             fclose(fid);
%         end
%         fid = fopen(fullfile(progressDir,strcat(jobName,'_percentComplete.txt')),'w');
%         fprintf(fid,'%f\n',percent);
%         fclose(fid);
%     end
    
    
    display(strcat('progress: ',num2str(numFinished),'/',num2str(nfiles)))
    filename = imagefiles{i};
    [pathstr,name,ext] = fileparts(filename);
    fnameOnly = strcat(name,ext);
    if exist(fullfile(output_dir,'failed',fnameOnly),'file')
        display('image already failed before, skipping')
        continue;
    end
    try
        I_Q = imread(filename);
    catch
        display('Image was in CMYK, skipping');
        failedImgs = failedImgs+1;
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
        continue
    end
    
    
    
    switch frontalization,
        case 'Hassner',
            [fidu_XY,frontal_raw,frontal_sym,hardsym_images]=demo(I_Q,detector,filename,[],faceSize,Model3D,'','','',eyemask);
            if ~isempty(fidu_XY) && ~isempty(frontal_raw)
                display ('Found face in image!')
            %         figure; imshow(I_Q); hold on; plot(fidu_XY(:,1),fidu_XY(:,2),'.');
            %         hold off; title('Query photo with detections overlaid');
            %         figure; imshow(frontal_raw); title('Frontalized no symmetry');
            %         figure; imshow(frontal_sym); title('Frontalized with soft symmetry');
            landmarks = fidu_XY;
            
            [full_dir,fname,ext] = fileparts(filename);
            base_dir = directoryPath;
            remaining_dir = full_dir(length(base_dir)+1:length(full_dir));
            display(remaining_dir);
            save_dir_sym = fullfile(output_dir,'sym',remaining_dir,[fname ext]);
            save_dir_asym = fullfile(output_dir,'asym',remaining_dir,[fname ext]);
            save_dir_failed = fullfile(output_dir,'failed',remaining_dir,[fname ext]);
           
            [full_dir,fname,ext] = fileparts(save_dir_sym);
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
            imwrite(frontal_sym,save_dir_sym);
            display(filename);
            display(save_dir_asym);
            else
                display('Could not detect face in image');
                failedImgs = failedImgs+1
                 try
                     display(save_dir_failed)
                imwrite(I_Q,save_dir_failed);
                catch
                end
            end
            %         display(fidu_XY)
        case 'Vito',
%             addpath('Vito') ;
            img = imresize(I_Q,1.1);
            frontal_raw = [];
%             try
            frontal_raw = viewmorph(img, models,dmodel, 'imsize',faceSize(1),'crop',2,'verbose',0,'detmethod',detector,'dlibdetector',dlibdetector,'dlibpredictor',dlibpredictor,'dlibmodule',dlibmod);
%             catch
%             end
            
            %             rawsize = length(frontal_raw)
%             display('made it here')
            if ~isempty(frontal_raw)
                display('frontalization succeeded!')
        
                 [full_dir,fname,ext] = fileparts(filename);
                base_dir = directoryPath;
                remaining_dir = full_dir(length(base_dir)+1:length(full_dir));
                save_dir_asym = fullfile(output_dir,'asym',remaining_dir,[fname ext]);
                save_dir_failed = fullfile(output_dir,'failed',remaining_dir,[fname ext]);
                [full_dir,fname,ext] = fileparts(save_dir_asym);
                if exist(full_dir, 'dir') == 0
                    mkdir(full_dir);
                end
                            [full_dir,fname,ext] = fileparts(save_dir_failed);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
                imwrite(frontal_raw,save_dir_asym);
                successImgs = successImgs + 1
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
                display('Could not detect face in image');
                failedImgs = failedImgs+1
                 try
                imwrite(I_Q,save_dir_failed);
                 catch
                 end
            end
        otherwise,
            disp('A combination of this Landmarker and Frontalization technique does not exist !!') ;
            
    end
    
    [full_dir,fname,ext] = fileparts(filename);
    base_dir = directoryPath;
    remaining_dir = full_dir(length(base_dir)+1:length(full_dir));
    save_dir_sym = fullfile(output_dir,'sym',remaining_dir,[fname ext]);
    save_dir_asym = fullfile(output_dir,'asym',remaining_dir,[fname ext]);
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
% parfor_progress(0,progressDir);
totalTime = toc;
disp('Finished all files!')
disp('Total time:')
disp(totalTime)
disp('Total images frontalized:')
disp(successImgs)
disp('Total images failed to frontalize:')
disp(failedImgs)

