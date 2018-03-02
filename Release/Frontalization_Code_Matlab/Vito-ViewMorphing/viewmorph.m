function im3 = viewmorph(im1,models,dmod,varargin)

% Required inputs
p = inputParser;
addRequired(p,'image',@isnumeric) % rgb image - matrix with the size of <nrow�ncol�3>
addRequired(p,'models',@isstruct)  %model structure (run function generateCMRmodel.m)

% Optional name-value pairs
addParameter(p,'verbose',1,@(x)x==0||x==1) %1/0 do/don't display figures
addParameter(p,'detmethod','CMR',@(x)any(strcmp(x,{'CMR','DLIB','ZR','INTRA'}))) %1/0- do/don't use intraface landmark detector to localize internal 49 facial landmarks
addParameter(p,'crop',1,@(x)x==1||x==2) %1-rectangular crop, 2-crop convex area bounded by landmarks
addParameter(p,'rollnorm',0,@(x)x==0||x==1) %1/0 - do/don't normalize roll angle (in-plane rotation)
addParameter(p,'prewarp',1,@(x)x==0||x==1) %1/0 do/don't perform prewarping
addParameter(p,'cpoints',0,@(x)x==0||x==1) %0 - automatically chosen control points, 1 - manually choose control points
addParameter(p,'transfun','new',@(x)any(strcmp(x,{'new','old'}))) %'old'/'new' - switch between old and new matlab built-in functions for spatial transformtations
addParameter(p,'postproc',2,@(x)x==0||x==1||x==2) %0 - no postporcessing, 1 - map morphed image to the average face, 2 - rectangular crop with fixed eye and mouth locations
addParameter(p,'imsize',[],@isnumeric) %size of the output image. If empty, don't resize.
addOptional(p,'dlibdetector',[])
addOptional(p,'dlibpredictor',[])
addOptional(p,'dlibmodule',[])
% Parse inputs
parse(p, im1, models, varargin{:});
im1 = p.Results.image;
model = p.Results.models;
verbose = p.Results.verbose;
detmethod = p.Results.detmethod;
crops = p.Results.crop;
rolln = p.Results.rollnorm;
prewarp = p.Results.prewarp;
cpoints = p.Results.cpoints;
mfun = p.Results.transfun;
postproc = p.Results.postproc;
imsize = p.Results.imsize;
ddetector = p.Results.dlibdetector;
dpredictor = p.Results.dlibpredictor;
dlibmodule = p.Results.dlibmodule;
%% Face detection and landmark localization
try
[lm1, bbox1] = detect_face_and_landmarks(im1,model,dmod,[],'detector',detmethod,'dlibdetector',ddetector,'dlibpredictor',dpredictor,'dlibmodule',dlibmodule);
catch
    im3 = [];
    return
end
if length(bbox1) > 0
if rolln %in-plane rotation to the upward position
    xl = mean(lm1(37:42,1));
    yl = mean(lm1(37:42,2));
    xr = mean(lm1(43:48,1));
    yr = mean(lm1(43:48,2));
    rollangle = atand((yr - yl)/(xr - xl));
    if 0 % rotate image and run both face detector and landmark localizator
        im1 = imrotate(im1,rollangle);
        display('redo')
        [lm1, bbox1] = detect_face_and_landmarks(im1,model,dmod,[],detmethod);
    else % recompute landmarks based on a roll angle
        theta = deg2rad(-rollangle);
        rotmat = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
        tform = affine2d(rotmat);
        [im1, lm1(:,2), lm1(:,1)] = transformImage(im1, lm1(:,2), lm1(:,1), tform);
        if 0 % run landmark detector on a rotated image for better localization
            bbox1 = [min(lm1(:,1))-(max(lm1(:,1))-min(lm1(:,1)))*.15,...
                     min(lm1(:,2))-(max(lm1(:,2))-min(lm1(:,2)))*.4,...
                     (max(lm1)-min(lm1))*1.3];
            lm1 = detect_face_and_landmarks(im1,model,bbox1,[],detmethod);
        end
    end    
end

im2 = fliplr(im1);
bbox2 = [size(im1,2)-bbox1(1)-bbox1(3) bbox1(2) bbox1(3:4)];
if 0 %run face and landmark detector on the flipped image
    [lm2, bbox2] = detect_face_and_landmarks(im2,model,dmod,bbox2,detmethod);
else %just recompute landmarks from the non-flipped image
    lm2 = flip_landmarks(lm1,size(im1));
end

switch crops
    
    case 1 %rectangular crop
        leftc = min(lm1(:,1));
        rightc = max(lm1(:,1));
        upc = min(lm1(:,2));
        downc = max(lm1(:,2));
        upc_ = round(upc - 0.9*(downc-upc));
        leftc_ = round(leftc - 0.1*(rightc-leftc));
        rightc_ = round(rightc + 0.1*(rightc-leftc));
        downc_ = round(downc + 0.1*(downc-upc));
        if upc_ < 1, upc_ = 1; end
        if leftc_ < 1, leftc_ = 1; end
        if downc_ > size(im1,1), downc_ = size(im1,1); end
        if rightc_ > size(im1,2), rightc_ = size(im1,2); end
        im1 = im1(upc_:downc_,leftc_:rightc_,:);   
        lm1 = lm1 - repmat([leftc_ upc_],size(lm1,1),1);
        bbox1(1,1) = bbox1(1,1) - leftc_;
        bbox1(1,2) = bbox1(1,2) - upc_;

        leftc = min(lm2(:,1));
        rightc = max(lm2(:,1));
        upc = min(lm2(:,2));
        downc = max(lm2(:,2));
        upc_ = round(upc - 0.9*(downc-upc));
        leftc_ = round(leftc - 0.1*(rightc-leftc));
        rightc_ = round(rightc + 0.1*(rightc-leftc)); 
        downc_ = round(downc + 0.1*(downc-upc));
        if upc_ < 1, upc_ = 1; end
        if leftc_ < 1, leftc_ = 1; end
        if downc_ > size(im2,1), downc_ = size(im2,1); end
        if rightc_ > size(im2,2), rightc_ = size(im2,2); end
        im2 = im2(upc_:downc_,leftc_:rightc_,:);   
        lm2 = lm2 - repmat([leftc_ upc_],size(lm2,1),1);
        bbox2(1,1) = bbox2(1,1) - leftc_;
        bbox2(1,2) = bbox2(1,2) - upc_;
        
    case 2 %convex crop
        im1_bw = zeros(size(im1,1),size(im1,2));
        im1_bw(sub2ind(size(im1_bw),lm1(:,2),lm1(:,1))) = 1;
        im1_bw = im1_bw == 1;
        im1_bw = bwconvhull(im1_bw);
        im1(~repmat(im1_bw,1,1,3)) = 0;
        im1( all(~im1_bw,2), :, :) = []; im1_bw(all(~im1_bw,2), :) = [];
        im1(:, all(~im1_bw,1), :) = [];
        bbox1(1:2) = double(bbox1(1:2)) - double(min(lm1));
        lm1 = lm1 - repmat(min(lm1),size(lm1,1),1) + 1;   
                
        im2_bw = zeros(size(im2,1),size(im2,2));
        im2_bw(sub2ind(size(im2_bw),lm2(:,2),lm2(:,1))) = 1;
        im2_bw = im2_bw == 1;
        im2_bw = bwconvhull(im2_bw);
        im2(~repmat(im2_bw,1,1,3)) = 0;
        im2( all(~im2_bw,2), :, :) = []; im2_bw(all(~im2_bw,2), :) = [];
        im2(:, all(~im2_bw,1), :) = [];
        bbox2(1:2) = double(bbox2(1:2)) - double(min(lm2));
        lm2 = lm2 - repmat(min(lm2),size(lm2,1),1) + 1;   
end

if verbose
    close all
    figure
    set(gcf,'Name','Detected landmarks')
    subplot 121
    imshow(im1)
    hold on
    rectangle('Position',bbox1(1,:),'EdgeColor','r')
    plot(lm1(:,1),lm1(:,2),'g.')
    text(lm1(:,1),lm1(:,2), cellstr(num2str(colon(1,length(lm1))')),'Color','g')
    subplot 122
    imshow(im2)
    hold on
    rectangle('Position',bbox2(1,:),'EdgeColor','r')
    plot(lm2(:,1),lm2(:,2),'g.')
    text(lm2(:,1),lm2(:,2), cellstr(num2str(colon(1,length(lm2))')),'Color','g')
end

% im1 = crop_eyes(im1,lm1);

%% Prewarp, morph and postwarp

lm1_ = lm1(18:end,:);
lm2_ = lm2(18:end,:);

if prewarp
    
    % Prewarp
    F = findfundmat_jk(lm1_, lm2_);
    [H1, H2] = computeH1H2(F);

    if strcmp(mfun,'old')
        tform1 = maketform('projective',H1');
        tform2 = maketform('projective',H2');
    else
        tform1 = projective2d(H1');
        tform2 = projective2d(H2');
    end

    % Get control points
    switch cpoints
        case 0
            idx = [37,46,49,55];
            failedTransform = 0;
            if crops == 2 %convex crop
                try
                [pw1, lm1_pw(:,2), lm1_pw(:,1)] = transformImage(im1, lm1(:,2), lm1(:,1), tform1);
                [pw2, lm2_pw(:,2), lm2_pw(:,1)] = transformImage(im2, lm2(:,2), lm2(:,1), tform2);
                catch
                    failedTransform = 1
                end
            else %rectangular crop
                try
                lm1_pw = [lm1; [1 1]; [1 size(im1,1)-1]; [size(im1,2)-1 1]; [size(im1,2)-1 size(im1,1)-1]];
                lm2_pw = [lm2; [1 1]; [1 size(im2,1)-1]; [size(im2,2)-1 1]; [size(im2,2)-1 size(im2,1)-1]];
                [pw1, lm1_pw(:,2), lm1_pw(:,1)] = transformImage(im1, lm1_pw(:,2), lm1_pw(:,1), tform1);
                [pw2, lm2_pw(:,2), lm2_pw(:,1)] = transformImage(im2, lm2_pw(:,2), lm2_pw(:,1), tform2);
                catch
                    failedTransform = 1
                end
            end
            if failedTransform == 1
                im3 = []
                return
            end
            X1 = lm1(idx,:);
            X2 = lm2(idx,:);
            UPW1 = lm1_pw(idx,:);
            UPW2 = lm2_pw(idx,:);

            im1tmp = imresize(im1, [size(pw1,1), size(pw1,2)]);
            im2tmp = imresize(im2, [size(pw2,1), size(pw2,2)]);
            pw1(repmat(sum(pw1,3)==0,1,1,3)) = im1tmp(repmat(sum(pw1,3)==0,1,1,3));
            pw2(repmat(sum(pw2,3)==0,1,1,3)) = im2tmp(repmat(sum(pw2,3)==0,1,1,3));

            if verbose
                figure
                set(gcf,'Name','Control points')
                subplot 221
                imshow(im1)
                title('Image 1')
                hold on
                plot(X1(:,1),X1(:,2),'g.')
                text(X1(:,1),X1(:,2), cellstr(num2str(colon(1,length(X1))')),'Color','g')
                subplot 222
                imshow(im2)
                title('Image 2')
                hold on
                plot(X2(:,1),X2(:,2),'g.')
                text(X2(:,1),X2(:,2), cellstr(num2str(colon(1,length(X2))')),'Color','g')
                subplot 223
                imshow(pw1)
                title('Prewarped image 1')
                hold on
                plot(UPW1(:,1),UPW1(:,2),'g.')
                text(UPW1(:,1),UPW1(:,2), cellstr(num2str(colon(1,length(UPW1))')),'Color','g')
                subplot 224
                imshow(pw2)
                title('Prewarped image 2')
                hold on
                plot(UPW2(:,1),UPW2(:,2),'g.')
                text(UPW2(:,1),UPW2(:,2), cellstr(num2str(colon(1,length(UPW2))')),'Color','g')
            end

        case 1
            if strcmp(mfun,'old')
                pw1 = imtransform(im1, tform1, 'size',size(im1));
                pw2 = imtransform(im2, tform2, 'size',size(im2));
            else
                pw1 = imwarp(im1, tform1);
                pw2 = imwarp(im2, tform2);
            end

            figure
            % Select corresponding 4 control points in image 1
            imshow(im1)
            X1 = ginput(4);
            % Select corresponding 4 control points in image 2
            imshow(im2)
            X2 = ginput(4);
            % Select corresponding 4 control points in prewarped image 1
            imshow(pw1);
            UPW1 = ginput(4);
            % Select corresponding 4 control points in prewarped image 2
            imshow(pw2);
            UPW2 = ginput(4);
    end

    % Morphing
    points_avg = lm1_pw*.5 + lm2_pw*.5;
    [points_avg,idx] = remove_duplicate_points(points_avg);
    tri = delaunayTriangulation(points_avg);
    [im,pts] = morph1(pw1,pw2,lm1_pw(idx,:),lm2_pw(idx,:),tri,.5,.5);
    
    % Post-warping
    c = .5;
    X = (1 - c)*X1 + c*X2;
    U = (1 - c)*UPW1 + c*UPW2;
    if size(U,1) ==  4 && size(X,1) == 4 && size(U,2) ==  2 && size(X,2) == 2
        try
    if strcmp(mfun,'old')
        tform = maketform('projective',U,X);
    else
        tform = fitgeotrans(U,X,'projective');
    end
        catch
            im3 = [];
            return
        end
    [im3, pts(:,2), pts(:,1)] = transformImage(im, pts(:,2), pts(:,1), tform);
    else
        im3 = [];
        return
    end
else
    
    %% Morph images without pre- and post- warping
    
    if crops==2 %convex crop
        lm1_ = lm1;
        lm2_ = lm2;
    else %rectangular crop
        lm1_ = [lm1; [1 1]; [1 size(im1,1)-1]; [size(im1,2)-1 1]; [size(im1,2)-1 size(im1,1)-1]];
        lm2_ = [lm2; [1 1]; [1 size(im2,1)-1]; [size(im2,2)-1 1]; [size(im2,2)-1 size(im2,1)-1]];
    end
    points_avg = lm1_*.5 + lm2_*.5;
    [points_avg,idx] = remove_duplicate_points(points_avg);
    tri = delaunayTriangulation(points_avg);
    [im3,pts] = morph1(im1,im2,lm1_(idx,:),lm2_(idx,:),tri,.5,.5);
    
end

%% Postprocessing
if size(pts,1) > 50 && size(pts,2) == 2

if postproc == 1
    
    [lm3,~] = detect_face_and_landmarks(im3,model,[],detmethod);
       
    lma = model.lm_avg;
    if size(lm3,1) == 66
        lma([61,65],:) = [];
        sclfak = ( mean((mean(lm3(43:48,1)) - mean(lm3(37:42,1)))./(mean(lma(43:48,1)) - mean(lma(37:42,1)))) + ...; 
            mean((mean(lm3(49:66,2)) - mean(lm3(37:48,2)))./(mean(lma(49:66,2)) - mean(lma(37:48,2)))) )/2; 
    else
        sclfak = ( mean((mean(lm3(43:48,1)) - mean(lm3(37:42,1)))./(mean(lma(43:48,1)) - mean(lma(37:42,1)))) + ...; 
            mean((mean(lm3(49:68,2)) - mean(lm3(37:48,2)))./(mean(lma(49:68,2)) - mean(lma(37:48,2)))) )/2; 
    end
    lma = round(lma*sclfak);
    lma = lma - repmat(mean(lma(37:48,:))-mean(lm3(37:48,:)),length(lm3),1);
    
    offset = round(mean(lm3(28:end,:)) - mean(lma(28:end,:)));
    lma(28:end,:) = lm3(28:end,:) - repmat(offset,length(lm3(28:end,:)),1);
    lma = lma([1:20,25:end],:);
    lm3 = lm3([1:20,25:end],:);

    [lma,idx] = remove_duplicate_points(lma);
    tri = generate_mapping_data(lma, 0);
    im3 = uint8(do_texture_mapping_with_delaunay(lm3(idx,:), im3, tri, 0));
    
    % Resize image
    if ~isempty(imsize)
        if length(imsize)>1 %width and hight given as input arguments
            im3 = imresize(im3,[imsize(1) imsize(2)]);
        else %one size parameter given as input argument
            im3 = imresize(im3,[1.1*imsize imsize]);
        end
    end

elseif postproc == 2
    
    eyedist = sqrt(sum((mean(pts(37:42,:)) - mean(pts(43:48,:))).^2));
    eyemouthdist = mean(pts([49,55],2)) - mean(pts(37:48,2));
    up = round(mean(pts(37:48,2)) - eyemouthdist*.55);
    down = round(mean(pts([49,55],2)) + eyemouthdist*.55);
    left = round(mean(pts(37:42,1)) - eyedist*.55);
    right = round(mean(pts(43:48,1)) + eyedist*.55);
   
    if up < 1, up = 1; end
    if down > size(im3,1), down = size(im3,1); end
    if left < 1, left = 1; end
    if right > size(im3,2), right = size(im3,2); end
    
    im3 = im3(up:down,left:right,:);
    
    % Resize image
    if ~isempty(imsize)
        if length(imsize)>1 %width and hight given as input arguments
            im3 = imresize(im3,[imsize(1) imsize(2)]);
        else %one size parameter given as input argument
            im3 = imresize(im3,[1.2*imsize imsize]);
        end
    end
    
end

%% Plot morphed image

if verbose
    figure%('units','normalized','outerposition',[0 0 1 1])
    set(gcf,'Name','View morphing')
    subplot 121
    imshow(im1,[])
    title('Original image')
    subplot 122
    imshow(im3)
    title('Morphed image')
    drawnow
end
else
    im3 = [];
end
else
    im3 = [];
end
