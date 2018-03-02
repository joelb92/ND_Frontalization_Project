function bbox = detect_faces_fp(X,VJdetector,num)

if nargin<3
    num=1;
end

% Init
bbox = [];
nsize = VJdetector.down_sampling_size;% = 320; % we want to shorter side of the image to 320 pixels

% Change size of Image
[a,b,~] = size(X);
if a<=b
    scale_fak = nsize/a; 
else
    scale_fak = nsize/b;
end
X=imresize(X,[round(a*scale_fak), round(b*scale_fak)],'bilinear');

% run detector
X2=uint8(imresize(X,[round(a*scale_fak), round(b*scale_fak)],'bilinear'));
boxes = VJdetector.fdetector.detect(X2, 'ScaleFactor',  1.1, ...
                            'MinNeighbors', 1, ...
                            'MinSize',      [20, 20]);                        

if num == 1  
    % here we run the profile detector only if the frontal one gets
    % nothing
    if numel(boxes)>=1 
       bbox = zeros(numel(boxes),4);
       for k=1:numel(boxes)
           rect =  boxes{k};
           bbox(k,:) = rect/scale_fak;
       end 
%         end
    else % if no frontsl faces are found, look for profile faces
       boxes1 = VJdetector.pdetector.detect(X2, 'ScaleFactor',  1.1, ...
                                'MinNeighbors', 1, ...
                                'MinSize',      [20, 20]);
       if numel(boxes1)>=1 
           bbox = zeros(numel(boxes1),4);
           for k=1:numel(boxes1)
               rect =  boxes1{k};
               bbox(k,:) = rect/scale_fak;
           end  
       end
    end
else
   % here we always run both detectors, but postprocess areas that are
   % the same
   if numel(boxes)>=1 
       bbox = zeros(numel(boxes),4);
       for k=1:numel(boxes)
           rect =  boxes{k};
           bbox(k,:) = rect/scale_fak;
       end
   end

   boxes1 = VJdetector.pdetector.detect(X2, 'ScaleFactor',  1.1, ...
                                'MinNeighbors', 1, ...
                                'MinSize',      [20, 20]);
  bbox1=[];
  if numel(boxes1)>=1 
       bbox1 = zeros(numel(boxes1),4);
       for k=1:numel(boxes1)
           rect =  boxes1{k};
           bbox1(k,:) = rect/scale_fak;
       end  
   end

   bboxf=bbox;
   for i=1:size(bbox1,1)
       x11 = bbox1(i,1);
       x12 = bbox1(i,1)+bbox1(i,3);
       y11 = bbox1(i,2);
       y12 = bbox1(i,2)+bbox1(i,4);
       flag = 1;
      for j=1:size(bbox,1)
          x21 = bbox(j,1);
          x22 = bbox(j,1)+bbox(j,3);
          y21 = bbox(j,2);
          y22 = bbox(j,2)+bbox(j,4);
          x_overlap = max(0, min(x12,x22) - max(x11,x21));
          y_overlap = max(0, min(y12,y22) - max(y11,y21));
          overlapArea = x_overlap * y_overlap;
          overlap_ratio = overlapArea/(bbox(j,3)*+bbox(j,3));
          if overlap_ratio<0.5
             flag =0;
          end
      end
      if flag
        bboxf = [bboxf;bbox1(i,:)];
      end
   end

   if ~isempty(bboxf)
       bbox = bboxf;
   else
       bbox=[];
   end

end
