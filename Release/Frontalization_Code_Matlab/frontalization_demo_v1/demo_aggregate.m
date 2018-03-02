function [  ] = demo_aggregate( image_file_path, landmarker, frontalization,AutoChoose)
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
predictor_file = fullfile(pwd,'shape_predictor_68_face_landmarks.dat')
detector = py.dlib.get_frontal_face_detector();
predictor = py.dlib.shape_predictor(predictor_file);
mod = py.importlib.import_module('dlib_detect_script_optimized');
mod = py.reload(mod);
str = sprintf('%s',[landmarker frontalization]) ;
I_Q = imread(image_file_path);
% if ~strcmp(image_file_path,'croppedtmp.jpg')
%     borderpx = 50;
% % bbox = bbox(1,:);
% [angle bbox] = estimatePoseZR(I_Q);
% maxdim = max(bbox(3),bbox(4));
% center = [bbox(1) + bbox(3)/2, bbox(2)+bbox(4)/2];
% maxdim = maxdim+borderpx*2;
% bbox = [center(1)-maxdim/2 center(2)-maxdim/2 maxdim maxdim];
% I_Q = imcrop(I_Q, bbox);
% imwrite(I_Q,'croppedtmp.jpg');
% demo_aggregate('croppedtmp.jpg',landmarker,frontalization,AutoChoose);
% end
% bbox = detect_faces_fp(I_Q,model.VJdetector,1);

frontal_sym = [];
switch frontalization,
    case 'Hassner',
        [fidu_XY,frontal_raw,frontal_sym,hardsym_images]=demo(I_Q,landmarker,image_file_path,[],[200,200],[],detector,predictor,mod);
%         figure; imshow(I_Q); hold on; plot(fidu_XY(:,1),fidu_XY(:,2),'.');
        %         hold off; title('Query photo with detections overlaid');
        %         figure; imshow(frontal_raw); title('Frontalized no symmetry');
        %         figure; imshow(frontal_sym); title('Frontalized with soft symmetry');
        landmarks = fidu_XY;
%         display(fidu_XY)
    case 'Vito',
        addpath('Vito') ;
        [warped_surface, fidu_XY,landmarked_img,frontal_raw,hardsym_images]=demoVito(I_Q,landmarker,image_file_path,[200,200],[],detector,predictor,mod);
        if ~isempty(landmarked_img)
        X=landmarked_img.img;
        landmarks=landmarked_img.landmarks;
        bbox=landmarked_img.bbox;
        
%         figure
%         imshow(X,[])
%         hold on
% %         rectangle('Position',bbox,'EdgeColor','r','LineWidth',2);
%         plot(landmarks(:,1)',landmarks(:,2)','go',...
%             'MarkerSize',2,...
%             'MarkerEdgeColor','g',...
%             'MarkerFaceColor','g')
%         display(fidu_XY)
%         figure
%         surf(warped_surface.X, warped_surface.Y, warped_surface.Z, warped_surface.C, 'edgecolor', 'none', 'FaceColor','texturemap')
        else
            display('Could not detect face in image');
        end
    otherwise,
        disp('A combination of this Landmarker and Frontalization technique does not exist !!') ;
        return
end
% [PATHSTR,NAME,EXT] = fileparts(image_file_path);
% savePath = fullfile(PATHSTR,NAME);
% saveName = strcat(NAME,'_',str);
% mkdir(strcat(fullfile(savePath)));

if ~isempty(hardsym_images)
if AutoChoose
    figure
    [angle bbox fidu_XY] = estimatePoseZR(I_Q);
    display(angle)
    main_img = [];
    if angle >= 45
        main_img = hardsym_images{1};
    elseif angle < 45 && angle >= 15
        main_img = uint8(frontal_sym);
        if isempty(main_img)
            [fidu_XY,frontal_raw,frontal_sym,hardsym_images]=demo(I_Q,landmarker,image_file_path,fidu_XY);
            main_img = frontal_sym;
        end
    elseif angle < 15 && angle > -15
        main_img = frontal_raw;
    elseif angle > -45 && angle <= -15
        main_img = uint8(frontal_sym);
        if isempty(main_img)
            [fidu_XY,frontal_raw,frontal_sym,hardsym_images]=demo(I_Q,landmarker,image_file_path,fidu_XY);
            main_img = frontal_sym;
        end
    elseif angle <= -45
        main_img = hardsym_images{2};
    end
    
    subplot(1,2,1)
    imshow(uint8(I_Q),[])
    title('Original')
    
    subplot(1,2,2)
    imshow(uint8(main_img),[])
    title('Best Frontalization')
%     imwrite(main_img,strcat(fullfile(savePath,saveName),'_best.jpg'));

    
else
    figure
    subplot(1,3,1)
    imshow(uint8(I_Q),[])
    hold on
    plot(landmarks(:,1)',landmarks(:,2)','go',...
            'MarkerSize',2,...
            'MarkerEdgeColor','g',...
            'MarkerFaceColor','g')
    title('Original')
    xlabel('[1]')
    hold off
    
    subplot(2,5,3)
    imshow(uint8(frontal_raw),[])
    title('No Symmetry')
    xlabel('[2]')
    
    subplot(2,5,4)
    if isempty(frontal_sym)
        frontal_sym = zeros(size(frontal_raw,1),size(frontal_raw,2));
    end
    imshow(uint8(frontal_sym),[])
    title('Soft Symmetry (if available)')
    xlabel('[3]')
    
    subplot(2,5,5)
    imshow(uint8(hardsym_images{1}),[])
    title('Hard Symmetry Right Side of Image [2]')
    if strcmp(frontalization,'Vito')
            title('Hard Symmetry Right Side')
    end
    xlabel('[4]')
    
    subplot(2,5,8)
    imshow(uint8(hardsym_images{2}),[])
    title('Hard Symmetry Left Side of Image [2]')
    if strcmp(frontalization,'Vito')
            title('Hard Symmetry Left Side')
    end
    xlabel('[5]')
    
    subplot(2,5,9)
    imshow(uint8(hardsym_images{3}),[])
    title('Hard Symmetry Right Side of Image [3]')
    if strcmp(frontalization,'Vito')
            title('Image [4] With Additional Center')
    end
%     xlabel('[6]')
%     subplot(2,5,10)
%     imshow(uint8(hardsym_images{4}),[])
%     title('Hard Symmetry Left Side of Image [3]')
%     if strcmp(frontalization,'Vito')
%             title('Image [5] With Additional Center')
%     end
    xlabel('[7]')
    
    % if length(hardsym_images) >4
    %     figure(7)
    % imshow(uint8(hardsym_images{5}),[])
    % title('With Hard Symmetry 2 candidate 2')
    % figure (8)
    % imshow(uint8(hardsym_images{6}),[])
    % title('With Hard Symmetry 2 candidate 2')
end
end
% else
%     imshow(I_Q)
%     title('No face was detected in image!')
end

