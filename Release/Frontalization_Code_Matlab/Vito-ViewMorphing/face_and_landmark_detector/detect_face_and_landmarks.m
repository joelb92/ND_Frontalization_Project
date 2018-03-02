function [lm, box] = detect_face_and_landmarks(im,model,dmod,varargin)

p = inputParser;
addRequired(p,'image',@isnumeric)
addRequired(p,'model',@isstruct)
addOptional(p,'dlibdetector',[])
addOptional(p,'dlibpredictor',[])
addOptional(p,'dlibmodule',[])
addOptional(p,'box',[],@isnumeric) % bounding box <x y w h>
addOptional(p,'detector','CMR',@(x)any(strcmp(x,{'CMR','DLIB','ZR','INTRA'}))) % detector method switch: 1-intraface, 2-our implementation

parse(p, im, model, varargin{:});
im = p.Results.image;
model = p.Results.model;
bbox = p.Results.box;
dmethod = p.Results.detector;
ddetector = p.Results.dlibdetector;
dpredictor = p.Results.dlibpredictor;
dmodule = p.Results.dlibmodule;
if strcmp(dmethod,'DLIB') || strcmp(dmethod,'ZR')
[model3D,lm1,box] = facial_feature_detection(dmethod,im,'',[],[],ddetector,dpredictor,dmodule,dmod);
if length(lm1) ~= 68
    box = [];
    lm = [];
%     display('empty')
else
%     display('not empty')
    box = double(box);
if  strcmp(dmethod,'ZR')
    lm = ZRtoDLIB(lm1);
else
    lm = lm1;
end
end
else


%% Detect face
extnd = round(max(size(im))/4);
ime = padarray(im,[extnd extnd],'both');
if isempty(bbox)
        box = detect_faces_fp(ime, model.VJdetector, 1);
        if size(box,1)>1, box = box(box(:,3)==max(box(:,3)),:); end
%         box = model.DM{1}.fd_h.detect(ime,'MinNeighbors',1,'ScaleFactor',1.1,'MinSize',[50 50]);
%         if length(box)>1, box = box{cellfun(@(c) c(3), box) == max(cellfun(@(c) c(3), box))}; else box = box{1}; end
else
    box = varargin{1};
    box(1:2) = box(1:2) + extnd;
end

%% Localize landmarks
lm1 = CMRfind_facial_landmarks(ime, box, model.Ldetector, 1);
if length(lm1) == 136
lm1 = reshape(lm1, [68,2]);
if strcmp(dmethod,'INTRA')
    lm2 = xx_track_detect(model,ime,box,struct('face_score',.3, 'min_neighbors',2, 'min_face_image_ratio',.15, 'compute_pose',1));
    lm2 = double(lm2.pred);
    lm = [lm1(1:17,:);lm2(1:end,:)];
else
    lm = lm1;
end
lm = lm - extnd;
box(1:2) = box(1:2) - extnd;
lm(lm<1) = 1;
lm(lm(:,1)>size(im,2),1) = size(im,2);
lm(lm(:,2)>size(im,1),2) = size(im,1);
else
    lm = [];
    box = [];
end
end
