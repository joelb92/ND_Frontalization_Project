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
function Frontalize_Hassner(image_dir,image_type,output_dir,detector,numCores,faceSize,progressDir,jobName)
% close all;
% clear all;
addpath calib;
addpath parfor_progress;
addpath kakearney-subdir/subdir;
% clusterName = parallel.importProfile(cluster_profile);

directoryPath = fullfile(image_dir);%'/Users/joel/Downloads/icbrw_Data/icbrw_ProbeImages/';
reducedFile = strcat('fileNamesReduced_',jobName,'.mat');
if ~exist(strcat('fileNames_',jobName,'.mat'), 'file')
    disp('loading all filenames...')
    folderLevel1 = dir(fullfile(directoryPath));
    imagefiles = [];
    
    for n = 1:length(folderLevel1)
        folder = folderLevel1(n);
        if folder.isdir
            imagefilestmp = dir(fullfile(directoryPath,folder.name,['*' image_type]));
            for i = 1:length(imagefilestmp)
                imagefilestmp(i).name = fullfile(directoryPath,folder.name,imagefilestmp(i).name);
            end
            imagefiles = [imagefiles ; imagefilestmp];
        end
    end
    
    save(strcat('fileNames_',jobName),'imagefiles');
    
else
    
    if exist(reducedFile,'file') == 2
        disp('loading filenames from reduced save file...')
        load(reducedFile);
    else
        disp('loading filenames from save file...')
        load(strcat('fileNames_',jobName));
    end
end

progressFileRegex = strcat('progress_',jobName,'*.txt');
progressFiles = dir(progressFileRegex);
alreadyProcessedIndexes = [];
display('Modifying array to reflect progress from progress files...')
for i = 1:length(progressFiles)
    display(['found progress file ' progressFiles(i).name])
    A = regexp( fileread(progressFiles(i).name), '\n', 'split');
    if length(A) >= 2
        start = str2num(A{2});
        stop = str2num(A{1});
        display(['Start: ' num2str(start) ' Stop: ' num2str(stop)])
        alreadyProcessedIndexes = [alreadyProcessedIndexes [start:stop]];
    end
end
if ~isempty(alreadyProcessedIndexes)
    display('Saving reduced file list')
    oldLength = length(imagefiles)
    imagefiles = imagefiles(setdiff(1:length(imagefiles),alreadyProcessedIndexes));
    newLength = length(imagefiles)
    save(reducedFile,'imagefiles');
    
end
delete(progressFileRegex);
disp('fould all files...');
% imagefiles = subdir(fullfile(directoryPath,['*' image_type]));
disp(['done reading in ' int2str(size(imagefiles)) ' files!'])
nfiles = length(imagefiles);    % Number of files found

% parpool('cvrl12',numCores)
% parpool('cvrl12',numCores)

poolobj = gcp('nocreate');

if length(poolobj) == 0  
if numCores > 0
    parpool('local',numCores)
else
    parpool('local')
end
else
    display('A local machine ')
end
% Load query image
eyemask = load('eyemask');
eyemask = eyemask.eyemask;
% load eyemask eyemask % mask to exclude eyes from symmetry
% load DataAlign2LFWa REFSZ REFTFORM % similarity transf. from rendered view to LFW-a coordinates

others = load('DataAlign2LFWa');
REFSZ = others.REFSZ;
REFTFORM = others.REFTFORM;

poolobj = gcp('nocreate');
Model3D = load(strcat('model3D',detector,'.mat'));
Model3D = Model3D.Model3D;
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
parfor_progress(nfiles,progressDir);
parfor i = 1:nfiles
    percent = parfor_progress(-1,progressDir);
    numFinished = numFinished+1;
    if 1 %keep track of progress
        t = getCurrentTask();
        progressFileName = strcat('progress_',jobName,'_',int2str(t.ID),'.txt');
        if exist(progressFileName,'file') == 2
            A = regexp( fileread(progressFileName), '\n', 'split');
            A{2} = sprintf('\n%i',i);
            fid = fopen(progressFileName, 'w');
            fprintf(fid, '%s', A{:});
            fclose(fid);
        else
            fid = fopen(progressFileName, 'w');
            fprintf(fid, '%i\n',i);
            fclose(fid);
        end
        fid = fopen(fullfile(progressDir,strcat(jobName,'_percentComplete.txt')),'w');
        fprintf(fid,'%f\n',percent);
        fclose(fid);
    end
    
    
    
    filename = imagefiles(i).name
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
    
    % I_Q = imresize(I_Q,.3);
    % load some data
    
    
    % Detect facial features with prefered facial feature detector
    
    % Note that the results in the paper were produced using SDM. We have found
    % other detectors to produce inferior frontalization results.
    if isempty(I_Q)
        display('could not load image')
        failedImgs = failedImgs+1;
    end
    [fidu_XY] = facial_feature_detection(detector,I_Q,filename);
    
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
    
    
    if isempty(fidu_XY)
        display('Failed to detect facial features / find face in image.');
        failedImgs = failedImgs+1;
        %     save_dir_sym = fullfile(output_dir,'sym_failed',remaining_dir,[fname ext]);
        %     save_dir_asym = fullfile(output_dir,'asym_failed',remaining_dir,[fname ext]);
        [full_dir,fname,ext] = fileparts(save_dir_failed);
        if exist(full_dir, 'dir') == 0
            mkdir(full_dir);
        end
        %     imwrite(I_Q,save_dir_sym);
        %     imwrite(I_Q,save_dir_asym);
        imwrite(I_Q,save_dir_failed);
    else
%         xcords = fidu_XY(:,1);
%         ycords = fidu_XY(:,2);
%         hull = convhull(xcords,ycords);
%         xcords = xcords(hull);
%         ycords = ycords(hull);
%         mask = repmat(uint8(poly2mask(xcords,ycords,size(I_Q,1),size(I_Q,2))),[1,1,3]);
%         I_Q = I_Q.*mask;
        % Estimate projection matrix C_Q
        if length(Model3D.threedee) == length(fidu_XY)
            successImgs = successImgs+1;
            [C_Q, est_A,~,~] = estimateCamera(Model3D, fidu_XY);
            threedee_model = Model3D;
            I_Q_3D = [I_Q()];
            % Render frontal view
            [frontal_sym, frontal_raw] = Frontalize(C_Q, I_Q, threedee_model.refU, eyemask);
            
            
            % Apply similarity transform to LFW-a coordinate system, for compatability
            % with existing methods and results
            frontal_sym = imtransform(frontal_sym,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);
            frontal_raw = imtransform(frontal_raw,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);
            frontal_sym = imresize(imcrop(frontal_sym.*mask,maskBB),faceSize);
            frontal_raw = imresize(imcrop(frontal_raw.*mask,maskBB),faceSize);
            
            try
                imwrite(frontal_sym,save_dir_sym);
                imwrite(frontal_raw,save_dir_asym);
                % disp(['wrote file: ', save_dir_sym]);
                % disp(['wrote file: ', save_dir_asym]);
            catch
                display('Could not write files, please check permissions!');
                failedImgs = failedImgs+1;
            end
        else
            display('Points in Model3D and detection are not equal')
            failedImgs = failedImgs+1;
            [full_dir,fname,ext] = fileparts(save_dir_failed);
            if exist(full_dir, 'dir') == 0
                mkdir(full_dir);
            end
            %     imwrite(I_Q,save_dir_sym);
            %     imwrite(I_Q,save_dir_asym);
            imwrite(I_Q,save_dir_failed);
        end
        % Display results
        
    end
    
end
parfor_progress(0,progressDir);
totalTime = toc;
disp('Finished all files!')
disp('Total time:')
disp(totalTime)
disp('Total images frontalized:')
disp(successImgs)
disp('Total images failed to frontalize:')
disp(failedImgs)

