% First frontalization tests - Vito


%% Init stuff
function [warped_surface, fid_XY, landmarked_img, frontal_raw,hardsym_images]=demoVito(X,detector,filename,faceSize,Model3D,dlibdetector,dlibpredictor,dlibmodule)
% load model
%clear all
%addpath(genpath('/Users/abharati/Desktop/HassnerS-Frontalization/frontalize.0.1.2/Vito/mexopencv-2.4')) % mexopencv - add your path
%addpath(genpath('C:\Users\VitoS\Dropbox\Face Alignment\demoForJoel'));
warped_surface = [];
fid_XY = [];
landmarked_img = [];
frontal_raw = [];
hardsym_images = {};
mpath = 'model/';
model = generateCMR_model1(mpath);

%estimate the pose of the person (to chose which side we should use for
%symmetry)


% some pictures
% pic1 = '/Users/abharati/Desktop/HassnerS-Frontalization/frontalize.0.1.2/Vito/pics/download.jpg';
% pic1 = '/Users/abharati/Desktop/HassnerS-Frontalization/frontalize.0.1.2/Vito/pics/images.jpg';
% pic1 = '/Users/abharati/Desktop/HassnerS-Frontalization/frontalize.0.1.2/Vito/pics/round-black.jpg';
% test_pic = pic1;

%% Do detection and landmarking and show results

% load picture
% X=imread(test_pic);

% do the detection
% addpath('/dlib/dlib/matlab')
% [bbox NA] = dlibFaceDetect(rgb2gray(X),[],1);
[model3D,landmarks,bbox] = facial_feature_detection(detector,X,filename,[],[],'','','');
if strcmp(detector,'ZhuRamanan')
    landmarks = ZRtoDLIB(landmarks);
end
% bbox = detect_faces_fp(X,model.VJdetector,1);
if length(bbox) > 0
    [C_Q, ~,~,~] = estimateCamera(model3D, landmarks);
    
  generalPose = SimplePoseDetector(X,landmarks)
% bbox=bbox(1,:);
% landmarks
% if strcmp(detector,'CMR')
%     % find the ladnmarks
%     landmarks = CMRfind_facial_landmarks(X,bbox,model.Ldetector,1);
%     landmarks=reshape(landmarks,[68,2]);
% elseif strcmp(detector,'dlib')
%         commandStr = ['python dlib_detect_script.py ',filename];
%         [status, commandOut] = system(commandStr);
%         if status==1
%              fprintf(commandOut);                                                                 
%         end
%         landmarks=[];
%         fid = fopen('landmarks.txt','r');
%         tline = fgetl(fid);
%         while ischar(tline)
%             strarr=strsplit(tline);
%             arr=[str2num(strarr{1,1}) str2num(strarr{1,2})];
%             landmarks=[landmarks;arr];
%             tline = fgetl(fid);
%         end 
%         fclose(fid);
%         landmarks = reshape(landmarks,68,2);
%         landmarks = double(landmarks);
% elseif strcmp(detector,'ZhuRamanan')
%         % Zhu and Ramanan detector [3]
%         % detect facial features on query
%         addpath(genpath('ZhuRamanan'))
%         load('ZhuRamanan/face_p146_small.mat','model');
%         model.interval = 5;
%         model.thresh = min(-0.65, model.thresh);
%         if length(model.components)==13 
%             posemap = 90:-15:-90;
%         elseif length(model.components)==18
%             posemap = [90:-15:15 0 0 0 0 0 0 -15:-15:-90];
%         else
%             error('Can not recognize this model');
%         end
%         
%         
%         I_Q_bs = detect(X, model, model.thresh);
%         if isempty(I_Q_bs)
%             return
%         end
%         I_Q_bs = clipboxes(X, I_Q_bs);
%         I_Q_bs = nms_face(I_Q_bs,0.3);
% 
%         if (isempty(I_Q_bs))
%             return;
%         end
%         x1 = I_Q_bs(1).xy(:,1);
%         y1 = I_Q_bs(1).xy(:,2);
%         x2 = I_Q_bs(1).xy(:,3);
%         y2 = I_Q_bs(1).xy(:,4);
%         fidu_XY = [(x1+x2)/2,(y1+y2)/2];
%         landmarks = double(fidu_XY);
% %         zhurama_hard_symmetry(I_Q,fidu_XY,Model3D);
% end
fid_XY = landmarks;
% rescale image and landmarks
% bbox_width = bbox(3);
% bbox_norm = 300;
% scale_fac = bbox_norm/bbox_width;
% bbox=bbox*scale_fac;
% landmarks=landmarks*scale_fac;
% X=imresize(X,scale_fac,'bilinear');

% figure
% imshow(X,[])
% hold on
% rectangle('Position',bbox,'EdgeColor','r','LineWidth',2);
% plot(landmarks(:,1)',landmarks(:,2)','go',...
%         'MarkerSize',2,...
%         'MarkerEdgeColor','g',...
%         'MarkerFaceColor','g')
landmarked_img.bbox=bbox;
landmarked_img.landmarks=landmarks;
landmarked_img.img=X;

%% Load the 3D model
sfname = '04676d119proc1.mat';
load([mpath,sfname]);
img = normalize8(img);

% this was used to annotate the model
% imshow(img,[])
% [x,y]=getpts();

pts = zeros(68,2);
pts(37,:) = [28.1718   79.1260];
pts(40,:) = [58.0954   78.5153];
pts(38,:) = [37.6374   74.5458];
pts(39,:) = [47.1031   73.0191];
pts(41,:) = [48.9351   84.6221];
pts(42,:) = [37.6374   84.6221];

pts(43,:) = [95.3473   77.5992];
pts(44,:) = [105.1183   73.9351];
pts(45,:) = [114.2786   74.2405];
pts(46,:) = [123.4389   79.7366];
pts(47,:) = [114.2786   82.7901];
pts(48,:) = [104.2023   82.7901];

pts(18,:) = [19.9275   68.4389];
pts(19,:) = [28.4771   61.4160];
pts(20,:) = [37.3321   59.2786];
pts(21,:) = [46.7977   61.7214];
pts(22,:) = [58.0954   65.3855];

pts(23,:) = [95.6527   66.9122];
pts(24,:) = [103.5916   61.4160];
pts(25,:) = [113.0573   58.3626];
pts(26,:) = [124.0496   60.8053];
pts(27,:) = [131.9885   66.3015];

pts(28,:) = [77.0267   73.0191];
pts(29,:) = [76.7214   88.8969];
pts(30,:) = [77.0267  102.6374];
pts(31,:) = [77.3321  117.2939];

pts(32,:) = [60.8435  125.8435];
pts(33,:) = [68.1718  127.6756];
pts(34,:) = [76.7214  129.2023];
pts(35,:) = [84.9656  126.7595];
pts(36,:) = [93.2099  125.2328];

pts(49,:) = [56.5687  153.0191];
pts(50,:) = [63.8969  148.1336];
pts(51,:) = [71.2252  145.9962];
pts(52,:) = [76.4160  148.7443];
pts(53,:) = [81.9122  145.6908];
pts(54,:) = [92.5992  148.7443];
pts(55,:) = [101.7595  155.4618];

pts(56,:) = [94.1260  159.1260];
pts(57,:) = [85.2710  164.0115];
pts(58,:) = [76.1107  164.3168];
pts(59,:) = [66.0344  163.0954];
pts(60,:) = [57.7901  159.1260];

pts(61,:) = [60.8435  153.0191];
pts(62,:) = [68.4771  152.4084];
pts(63,:) = [74.5840  153.6298];
pts(64,:) = [81.6069  152.4084];
pts(65,:) = [90.1565  153.9351];
pts(66,:) = [82.8282  157.2939];
pts(67,:) = [74.8893  156.9885];
pts(68,:) = [68.7824  155.4618];

pts(1,:) = [6.1870   83.7061];
pts(2,:) = [8.0191  103.8588];
pts(3,:) = [11.3779  122.4847];
pts(4,:) = [16.5687  140.5000];
pts(5,:) = [24.5076  156.9885];
pts(6,:) = [36.4160  171.9504];
pts(7,:) = [47.1031  184.1641];
pts(8,:) = [61.7595  193.3244];
pts(9,:) = [77.0267  198.2099];
pts(10,:) = [93.6527  194.5458];
pts(11,:) = [111.8359  184.1641];
pts(12,:) = [121.8282  172.8664];
pts(13,:) = [131.6832  159.1260];
pts(14,:) = [138.4008  142.0267];
pts(15,:) = [141.7595  120.9580];
pts(16,:) = [145.4237  100.1947];
pts(17,:) = [144.8130   80.3473];

% we mirror the 3D model and coordinates to make them symetric
eye_mirror_perms = [17:-1:1,27:-1:18,28:31,36:-1:32,46:-1:43,48,47,40:-1:37,42,41,55:-1:49,60:-1:56,65:-1:61,68:-1:66];
pts2 = [150-pts(:,1),pts(:,2)];
pts2 = pts2(eye_mirror_perms,:);

pts = round((pts+pts2)/2);

% We add a few additional coordinates to make the facial area bigger (not really needed)
ptsx = [pts];
pts1(1,:) = [22.0000   14.0000];
pts1(2,:) = [75.0000    1.0000];
pts1(3,:) = [128.0000   14.0000];
pts = [pts;pts1];

% produce the Z-values of the annotated points
for i=1:68
    Z(i)=img(round(pts(i,2)),round(pts(i,1)));
end


% assemble 3D points into homogenous form
x3d = zeros(68,3);
homogen3D = zeros(68,4);
for i=1:68
    x3d(i,:) =[round(pts(i,1)), round(pts(i,2)), Z(i)];
    homogen3D(i,:) = [round(pts(i,1)), round(pts(i,2)), Z(i), 1];
end

%normalize 3D points to origin at (0,0,0) and average distance from origin
%sqrt(3) - this is important for camera estimation
mua = mean(x3d);
x3d_zm = x3d-repmat(mua,size(x3d,1),1);
d_sq = sum((x3d_zm(:,1).^2),2);
mean_sq = mean(d_sq);
s = sqrt(mean_sq/3);
scalexyz = sqrt(  (sum(sum((x3d_zm.^2),2)))/(3*68)   );
scale = [1/scalexyz, 1/scalexyz,1/scalexyz];

U = [1/scalexyz,0,0;0,1/scalexyz,0;0 ,0,1/scalexyz;0,0,0];
U = [U,[-(mua.*scale)';1]];
homo_x3d = [x3d,ones(68,1)]';
normX3d = U*homo_x3d;

% assemble 2D points in homogeneuos form
x2d = zeros(68,2);
for i=1:68
    x2d(i,:) =[round(landmarks(i,1)), round(landmarks(i,2))];
end
x2d = x2d([8:10,18:68],:); % exclude face outline

%normalize 2D points to origin at (0,0) and average distance from origin
%sqrt(2)
mua = mean(x2d);
x2d_zm = x2d-repmat(mua,size(x2d,1),1);
d_sq = sum((x2d_zm(:,1).^2),2);
mean_sq = mean(d_sq);
s = sqrt(mean_sq/2);
scalexyz = sqrt(  (sum(sum((x2d_zm.^2),2)))/(2*size(x2d,1))   );
scale = [1/scalexyz, 1/scalexyz];
T = [1/scalexyz,0;0,1/scalexyz;0,0];
T = [T,[-(mua.*scale)';1]];
homo_x2d = [x2d,ones(size(x2d,1),1)]';
normX2d = T*homo_x2d;

% compute camera
zero_vec = zeros(1,4);
X_3D = zeros(136,8);
cont=1;
for i=1:2:136
    vec3D = normX3d(:,cont)';
    X_3D(i,:) = [vec3D zero_vec];
    X_3D(i+1,:) = [zero_vec vec3D];
    cont=cont+1;
end
X_3D = X_3D([15:20,35:136],:);


x_2D = normX2d(1:2,:);
x_2D = x_2D(:);

X_3Ds = X_3D;
x_2Ds = x_2D;
Ptv = (X_3Ds'*X_3Ds)\X_3Ds'*x_2Ds;

Pt = [reshape(Ptv,[4,2])';[0 0 0 1]];
P = T\Pt*U; % this is our final camera


% compute affine 3D to 2D camera
x_proj = P*homogen3D'; % here we see how good are comera is


homogen3Dall = zeros(71,4);
for i=1:71
    Z(i)=img(round(pts(i,2)),round(pts(i,1)));
    if isnan(Z(i))
        Z(i)=0;
    end
    homogen3Dall(i,:) = [round(pts(i,1)), round(pts(i,2)), Z(i), 1];
end




% compute triangulation and frontalize
% triangulation_pts = homogen3Dall(:,1:2);
% tri = generate_3Dtriangulation_mapping_data(triangulation_pts,img,0);
try
tri = generate_3Dtriangulation_mapping_data(ptsx,img,0);

% figure
% subplot(1,2,1)
% title('Original')
% imshow(uint8(X),[])
% subplot(1,2,2)
% title('3D Warpped')
% imshow(uint8(mapped_img),[])

mapped_img = do_3Dtexture_mapping_with_delaunay(landmarks, X, img, tri,P,0);
catch ME
    warped_surface = [];
    landmarked_img = [];
    frontal_raw = [];
    hardsym_images = [];
    return
end
mapped_img = imcrop(mapped_img,[1 51 size(mapped_img,2) size(mapped_img,1)-50]);

mapped_img1 = fliplr(mapped_img);
% frontal_raw = imcrop(uint8(mapped_img),[1 51 size(mapped_img,2) size(mapped_img,1)-50]);

frontal_raw = uint8(mapped_img);
maped1 = uint8([mapped_img(:,1:75,:), mapped_img1(:,76:end,:)]);
maped2 = uint8([mapped_img1(:,1:75,:), mapped_img(:,76:end,:)]);
maped3 = uint8([mapped_img(:,1:76,:), mapped_img1(:,74:end,:)]);
maped4 = uint8([mapped_img1(:,1:76,:), mapped_img(:,76:end,:)]);
maped5 = uint8([mapped_img(:,1:77,:), mapped_img1(:,73:end,:)]);
maped6 = uint8([mapped_img1(:,1:77,:), mapped_img(:,73:end,:)]);

maped1 = imresize(maped1,faceSize);
maped2 = imresize(maped2,faceSize);
maped3 = imresize(maped3,faceSize);
maped4 = imresize(maped4,faceSize);
maped5 = imresize(maped5,faceSize);
maped6 = imresize(maped6,faceSize);

if generalPose == 1
    hardsym_images = {maped1,maped3,maped5};
elseif generalPose == -1
    hardsym_images = {maped2,maped4,maped6};
else
    hardsym_images = {maped1,maped2,frontal_raw};
end



% figure
% subplot(2,4,1)
% title('Original')
% imshow(uint8(X),[])
% subplot(2,4,2)
% imshow(uint8(mapped_img),[])
% title('3D Warpped')
% 
% subplot(2,3,3)
% imshow(uint8(maped1),[])
% title('3D Warpped - hard symmetric 1')
% 
% subplot(2,3,4)
% imshow(uint8(maped2),[])
% title('3D Warpped - hard symmetric 1')
% 
% subplot(2,3,5)
% imshow(uint8(maped3),[])
% title('3D Warpped - hard symmetric 1')
% 
% subplot(2,3,6)
% imshow(uint8(maped4),[])
% title('3D Warpped - hard symmetric 1')



% subplot(2,4,7)
% maped5 = [mapped_img(:,1:77,:), mapped_img1(:,73:end,:)];
% imshow(uint8(maped5),[])
% title('3D Warpped - hard symmetric 1')
% 
% subplot(2,4,8)
% maped6 = [mapped_img1(:,1:77,:), mapped_img(:,73:end,:)];
% imshow(uint8(maped6),[])
% title('3D Warpped - hard symmetric 1')


% this produces the plot of the 3D template with the mapped texture

warped_surface.X=1:150;
warped_surface.Y=1:200;
warped_surface.Z=img;
warped_surface.C=uint8(mapped_img);
warped_surface.edgeColor='none';
warped_surface.FaceColor='texturemap';
frontal_raw = imresize(frontal_raw,faceSize);

else
    warped_surface = [];
    landmarked_img = [];
    frontal_raw = [];
    hardsym_images = [];
    return
end
%figure
%surf(1:150, 1:200, img, uint8(mapped_img), 'edgecolor', 'none', 'FaceColor','texturemap')


end


