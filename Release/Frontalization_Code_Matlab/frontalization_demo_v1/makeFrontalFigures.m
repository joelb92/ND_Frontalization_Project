function [  ] = makeFrontalFigures( image_file_path, landmarker, frontalization)
% This function takes in the location of the image, the techniques of
% landmarker, that can be one of the following :
% -'dlib'
% -'ZhuRamanan'
% -'CMR'
%
% the possible options of frontalization are :
% -'Hassner'
% -'Vito'
%   Detailed explanation goes here
close all
str = sprintf('%s',[landmarker frontalization]) ;
I_Q = imread(image_file_path);
% imSize = size(I_Q)
% scaleRatio = 1000.0/double(imSize(1))
% I_Q = imresize(I_Q,scaleRatio);
frontal_sym = [];
switch frontalization,
    case 'Hassner',
        landmarker='dlib';
        [fidu_XY,frontal_raw,frontal_sym,hardsym_images]=demo(I_Q,landmarker,image_file_path,[]);
        figure; imshow(I_Q); hold on; plot(fidu_XY(:,1),fidu_XY(:,2),'.');
%         hold off; title('Query photo with detections overlaid');
%         figure; imshow(frontal_raw); title('Frontalized no symmetry');
%         figure; imshow(frontal_sym); title('Frontalized with soft symmetry');
        
        
    case 'Vito',
        addpath('Vito') ;
        [warped_surface, landmarked_img,frontal_raw,hardsym_images]=demoVito(I_Q,landmarker,image_file_path);
        
%         X=landmarked_img.img;
%         landmarks=landmarked_img.landmarks;
%         bbox=landmarked_img.bbox;
%         
%         figure
%         imshow(X,[])
%         hold on
%         rectangle('Position',bbox,'EdgeColor','r','LineWidth',2);
%         plot(landmarks(:,1)',landmarks(:,2)','go',...
%             'MarkerSize',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','g')
%         
%         figure
%         surf(warped_surface.X, warped_surface.Y, warped_surface.Z, warped_surface.C, 'edgecolor', 'none', 'FaceColor','texturemap')
    otherwise,
        disp('A combination of this Landmarker and Frontalization technique does not exist !!') ;
        return
end
% figure
% subplot(2,4,1)
imshow(uint8(I_Q),[])
[PATHSTR,NAME,EXT] = fileparts(image_file_path);
savePath = fullfile(PATHSTR,NAME);
saveName = strcat(NAME,'_',str);
mkdir(strcat(fullfile(savePath)));

if ~isempty(frontal_raw)
title('Original')
figure(1)
% subplot(2,4,2)
% imshow(uint8(frontal_raw),[])
title('Without Symmetry')
% savefig(strcat(fullfile(savePath,saveName),'_raw'));
imwrite(uint8(frontal_raw),strcat(fullfile(savePath,saveName),'_raw.jpg'));
end

if ~isempty(frontal_sym)
% subplot(2,4,3)
figure(2)
% imshow(uint8(frontal_sym),[])
title('With Soft Symmetry (if available)')
% savefig(strcat(fullfile(savePath,saveName),'_soft'));
imwrite(frontal_sym,strcat(fullfile(savePath,saveName),'_soft.jpg'));
end

if ~isempty(frontal_raw)
% subplot(2,4,4)
figure(3)
% imshow(uint8(hardsym_images{1}),[])
title('With Hard Symmetry Side 1')
% savefig(strcat(fullfile(savePath,saveName),'_hard1'));
imwrite(uint8(hardsym_images{1}), strcat(fullfile(savePath,saveName),'_hard1.jpg'));

% subplot(2,4,5)
figure(4)
% imshow(uint8(hardsym_images{2}),[])
title('With Hard Symmetry Side 2')
% savefig(strcat(fullfile(savePath,saveName),'_hard2'));
imwrite(uint8(hardsym_images{2}),strcat(fullfile(savePath,saveName),'_hard2.jpg'));


figure(5)
% subplot(2,4,6)
% imshow(uint8(hardsym_images{3}),[])
title('With Hard Symmetry Side 1 cadidate 2')
% savefig(strcat(fullfile(savePath,saveName),'_hard3'));
imwrite(uint8(hardsym_images{3}),strcat(fullfile(savePath,saveName),'_hard3.jpg'));


figure(6)
% subplot(2,4,7)
% imshow(uint8(hardsym_images{4}),[])
title('With Hard Symmetry 2 candidate 2')
% savefig(strcat(fullfile(savePath,saveName),'_hard4'));
imwrite(uint8(hardsym_images{4}),strcat(fullfile(savePath,saveName),'_hard4.jpg'));


if length(hardsym_images) >4
    figure(7)
% subplot(2,4,7)
% imshow(uint8(hardsym_images{5}),[])
title('With Hard Symmetry 2 candidate 2')
% savefig(strcat(fullfile(savePath,saveName),'_hard5'));
imwrite(uint8(hardsym_images{5}),strcat(fullfile(savePath,saveName),'_hard5.jpg'));

figure (8)
% subplot(2,4,7)
% imshow(uint8(hardsym_images{6}),[])
title('With Hard Symmetry 2 candidate 2')
% savefig(strcat(fullfile(savePath,saveName),'_hard6'));
imwrite(uint8(hardsym_images{6}),strcat(fullfile(savePath,saveName),'_hard6.jpg'));
end
end

end

