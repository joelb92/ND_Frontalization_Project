function [side] = SimplePoseDetector(I_Q,landmarks)
   convinds = convhull(landmarks(:,1),landmarks(:,2));
   faceHull = [];
   for i = 1:length(convinds)
       faceHull = [faceHull; landmarks(convinds(i),:)];
   end
   faceMask = poly2mask(faceHull(:,1),faceHull(:,2),size(I_Q,1),size(I_Q,2));
   props = regionprops(faceMask,'BoundingBox');
   bbox = props.BoundingBox;
   cX = 0;
   cY = 0;
   centerPoints = 28:36;
   for i = 1:9
       cX = cX+landmarks(centerPoints(i),1);
       cY = cY+landmarks(centerPoints(i),2);
   end
   cX = cX/9;
   cY= cY/9;
   
   fcX = 0;
   fcY = 0;
   for i = 1:size(faceHull,1)
       fcX = fcX+faceHull(i,1);
       fcY = fcY+faceHull(i,2);
   end
   fcX = fcX/size(faceHull,1);
   fcY= fcY/size(faceHull,1);
   dist = cY-fcY;
   if(abs(cX-fcX) <= bbox(3)*.15) 
       side = 0;
   elseif (cX-fcX > 0) 
       side = 1;
   else
       side = -1;
   end
   return
%  ACC_CONST = 800; 
%     I_Q = double(I_Q);
% 
%     bgind = sum(abs(refU),3)==0;
% 
%     % count the number of times each pixel in the query is accessed
%     threedee = reshape(refU,[],3)';
%     tmp_proj = C_Q * [threedee;ones(1,size(threedee,2))];
%     tmp_proj2 = tmp_proj(1:2,:)./ repmat(tmp_proj(3,:),2,1);
%     
% 
%     bad = min(tmp_proj2)<1 | tmp_proj2(2,:)>size(I_Q,1) | tmp_proj2(1,:)>size(I_Q,2) | bgind(:)';
%     tmp_proj2(:,bad) = [];
% 
%     ind = sub2ind([size(I_Q,1),size(I_Q,2)], round(tmp_proj2(2,:)),round(tmp_proj2(1,:)));
% 
%     synth_frontal_acc = zeros(size(refU,1),size(refU,2));
%     
%     ind_frontal = 1:(size(refU,1)*size(refU,2));
%     ind_frontal(bad) = [];
%         
%     [c,~,ic] = unique(ind);
%     count = hist(ind,c);
%     synth_frontal_acc(ind_frontal) = count(ic);
% 
%     synth_frontal_acc(bgind) = 0;
%     synth_frontal_acc = imfilter(synth_frontal_acc,fspecial('gaussian', 16, 30),'same','replicate');
%     
%     % create synthetic view, without symmetry
%     
%     % which side has more occlusions?
%     midcolumn = round(size(refU,2)/2);
%     sumaccs = sum(synth_frontal_acc);
%     sum_left = sum(sumaccs(1:midcolumn));
%     sum_right = sum(sumaccs(midcolumn+1:end));
%     sum_diff = sum_left - sum_right;
%     
%     if abs(sum_diff)>ACC_CONST % one side is occluded
%         if sum_diff > ACC_CONST % left side of face has more occlusions
%             side =  1
%             return
%         else % right side of face has occlusions
%             side = -1  
%             return
%         end
%     else
%         side =  0
%         return
    end