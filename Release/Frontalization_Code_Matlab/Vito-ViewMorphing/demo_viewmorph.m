% Demonstration of the View-morphing method

% addpath('C:\opencv\mexopencv')
% addpath(genpath('C:\Users\VitoS\Documents\Administracija\Consulting\Jobs\EngineUpgrade_20141119\working_copy_code\opencv'))
% cdir = fileparts(mfilename('fullpath'));
% addpath([cdir,'\face_and_landmark_detector'])
addpath('face_and_landmark_detector/')
% Load models
models = generateCMRmodel('models/');

% Load images
imgs = dir('pics/probe/*.jpg');
for c = 1 : length(imgs)
    img = imresize(imread(['pics/probe/',imgs(c).name]),1.1);
    img_morph = viewmorph(img, models, 'imsize',200,'crop',2);
    imwrite(img_morph,['pics/processed/',imgs(c).name])
    disp('Press any button to continue ...')
    pause
end
