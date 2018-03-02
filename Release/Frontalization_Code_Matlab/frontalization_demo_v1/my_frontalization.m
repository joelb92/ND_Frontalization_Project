
clear;
close all;
addpath calib

% Load query image
images=dir('/Users/abharati/Desktop/reposing-images/JPEG/probe/*.jpg');

for ii=1:length(images)
    location='/Users/abharati/Desktop/reposing-images/JPEG/probe/';
    I_Q =imresize(imread(strcat(location,images(ii).name)),[250 250]);
    filename=images(ii).name;
    strarr=strsplit(filename,'.');
    %I_Q = imread('test.jpg');
    % load some data 
    dlibLoc='/Users/abharati/Desktop/1st Sem/dlib-18.17/examples/build1/JPEG/';
    filename=strcat(strarr{1,1},'0.txt');
    createDlibLandmarks(filename)
    load eyemask eyemask % mask to exclude eyes from symmetry
    load DataAlign2LFWa REFSZ REFTFORM % similarity transf. from rendered view to LFW-a coordinates

    % Detect facial features with prefered facial feature detector 
    detector = 'dlib'; % alternatively 'ZhuRamanan', 'dlib'
    % Note that the results in the paper were produced using SDM. We have found
    % other detectors to produce inferior frontalization results. 
    fidu_XY = [];
    pause(3)
    facial_feature_detection;
    if ~isempty(fidu_XY)
        

    % Estimate projection matrix C_Q
    [C_Q, ~,~,~] = estimateCamera(Model3D, fidu_XY);

    % Render frontal view
    [frontal_sym, frontal_raw] = Frontalize(C_Q, I_Q, Model3D.refU, eyemask);


    % Apply similarity transform to LFW-a coordinate system, for compatability
    % with existing methods and results
    frontal_sym = imtransform(frontal_sym,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);
    frontal_raw = imtransform(frontal_raw,REFTFORM,'XData',[1 REFSZ(2)], 'YData',[1 REFSZ(1)]);


    % Display results
    % figure; imshow(I_Q); title('Query photo');
    % figure; imshow(I_Q); hold on; plot(fidu_XY(:,1),fidu_XY(:,2),'.'); hold off; title('Query photo with detections overlaid');
    % figure; imshow(frontal_raw); title('Frontalilzed no symmetry');
    % figure; imshow(frontal_sym); title('Frontalilzed with soft symmetry');
    pause(3)
    saveLoc='/Users/abharati/Desktop/reposing-images/JPEG/dlib_frontalized_probe/';
    imwrite(frontal_sym, strcat(saveLoc,images(ii).name));
    else
        images(ii).name
    end

end
