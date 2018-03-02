function [lm, box] = detect_face_and_landmarks(im,model,varargin)
global box_buf
p = inputParser;
addRequired(p,'image',@isnumeric)
addRequired(p,'model',@isstruct)
addOptional(p,'box',[],@isnumeric) % Bounding box <x y w h>
addOptional(p,'fdetect',4,@(x)any(x==[0,1,2,3,4])) % Face detector method
addOptional(p,'ldetect',1,@(x)x==0||x==1) % Landmark detector method switch: 1-intraface, 0-our implementation

parse(p, im, model, varargin{:});
im = p.Results.image;
model = p.Results.model;
box = p.Results.box;
fdetect = p.Results.fdetect;
ldetect = p.Results.ldetect;

%% Detect face
extnd = round(max(size(im))/4);
ime = padarray(im,[extnd extnd],'both');
if isempty(box)
    switch fdetect %face detector switch
        case 1 % 1st option (Matlab Viola-Jones implementation)
            box = step(vision.CascadeObjectDetector('FrontalFaceCART'), ime);
            if isempty(box), box = step(vision.CascadeObjectDetector('FrontalFaceLBP'), ime); end
            if isempty(box), box = step(vision.CascadeObjectDetector('ProfileFace'), ime); end
            if size(box,1)>1, box = box(box(:,3)==max(box(:,3)),:); end
        case 2 % 2nd option
            box = detect_faces_fp(ime, model.VJdetector, 1);
            if size(box,1)>1, box = box(box(:,3)==max(box(:,3)),:); end
        case 3 % 3rd option
            box = model.DM{1}.fd_h.detect(ime,'MinNeighbors',1,'ScaleFactor',1.1,'MinSize',[25 25]);
            if length(box)>1, box = box{cellfun(@(c) c(3), box) == max(cellfun(@(c) c(3), box))}; else box = box{1}; end
        case 4 % 4th option (http://www.cbsr.ia.ac.cn/users/scliao/projects/npdface/)
            box = DetectFace(model.fd_frontal, ime);
            if isempty(box), box = DetectFace(model.fd_unconstrain, ime); end
            if length(box)>1, box = box([box.size] == max([box.size])); end
            if ~isempty(box)
                box = [box(1).col box(1).row box(1).size box(1).size];
                scl = 1.4;
                box = [box(1:2)-(1-scl/2)*box(3) box(3:4)*scl];
            end
            
    end 
    if isempty(box) %If face detector fails, use the central part of the image as a face region.
       warning('Couldn''t detect any faces, using the middle part of the image as a face detection box.')
       extnd = 0;
       ime = im;
       box = [size(im,1)/4 size(im,2)/4 min([size(im,1) size(im,2)])/2 min([size(im,1) size(im,2)])/2];
       %box = [size(im,1)/5 size(im,2)/5 size(im,1)*3/5 size(im,2)*3/5];
    else
        box_buf = [box_buf; box];
    end
else
    box = varargin{1};
    box(1:2) = box(1:2) + extnd;
end

%% Localize landmarks
lm1 = CMRfind_facial_landmarks(ime, box, model.Ldetector, 1);
lm1 = reshape(lm1(1,:), [68,2]);
if ldetect == 1 %Intraface
    lm2 = xx_track_detect(model,ime,box,struct('face_score',0, 'min_neighbors',2, 'min_face_image_ratio',.15, 'compute_pose',1));
    lm2 = double(lm2.pred);
    if size(lm2,1)==49
        lm = [lm1(1:17,:);lm2(1:end,:)];
        %lm = [lm1(1:17,:);lm2(1:43,:);lm1(61,:);lm2(44:46,:);lm1(65,:);lm2(47:end,:)];
    else
        lm = lm1;
    end
else
    lm = lm1;
end
lm = round(lm - extnd);
box(1:2) = box(1:2) - extnd;
lm(lm<1) = 1;
lm(lm(:,1)>size(im,2),1) = size(im,2);
lm(lm(:,2)>size(im,1),2) = size(im,1);
