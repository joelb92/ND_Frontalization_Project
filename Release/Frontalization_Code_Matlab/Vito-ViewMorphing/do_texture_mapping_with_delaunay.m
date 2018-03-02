% The function performs texture mapping based on delaunay triangualtion. 
% 
% Prototype:
%     mapped_img = 
%     do_texture_mapping_with_delaunay(probe_mesh_pts, prob_image, tri,verbose)
%     
% Inputs: 
%     probe_mesh_pts ... a Nx2 matrix of coordiantes that define the probe shape
%     prob_image ... the probe image to be warped
%     tri ... structure precomputed with "generate_mapping_data"
%     verbose ... report to prompt (binary)
%     
% Output:
%     mapped_img ... an image with warped image texture 
% 
% The function takes as the first input a set of reference mesh points that 
% define the target shape, to which we want to map our texture. The second 
% input is the set of mesh points in the probe image. If we assume that the 
% mesh points correspond (they are properly ordered) and that the same or a 
% similar object is present on the reference and probe images, then this 
% function can be used to create similar appearences of the probe and 
% reference images. 
% 
% The function first computes a delauney triangulation of the reference mesh 
% and then mapps each triangle in the probe image to the corresponding 
% triangle in the refernce shape.
% 
% Example:
%     An example of the use of this function would be to align facial images 
%     to a predefined target shape.
%     
% Author: Vitomir Struc
% Copyright, FE UL, LUKS, 2014
% Date: 29. 8. 2014

function mapped_img = do_texture_mapping_with_delaunay(probe_mesh_pts, prob_image, tri,verbose)

%% init
if nargin<4
    verbose = 0;
end

probe_mesh_pts=round(probe_mesh_pts);

%% clean probe image
max_x = max(probe_mesh_pts(:,1));
min_x = min(probe_mesh_pts(:,1));
max_y = max(probe_mesh_pts(:,2));
min_y = min(probe_mesh_pts(:,2));
size_y = max_y-min_y+1;
size_x = max_x-min_x+1;

[ax,bx,cx]=size(prob_image);
if cx==1
    prob_img_clean = prob_image(min_y:max_y,min_x:max_x);
    [mn nn] = size(prob_img_clean);
else
    prob_img_clean = prob_image(min_y:max_y,min_x:max_x,:);
    [mn nn dd] = size(prob_img_clean);
end

%% normalize probe shape coordinates and find triangle points
probe_pts = [probe_mesh_pts(:,1)-min_x+1,probe_mesh_pts(:,2)-min_y+1];
% [In,Jn] = ind2sub([mn nn],1:(mn*nn));
% 
% prob_inds = pointLocation(tri.dres,[Jn;In]');
% prob_inds = sparse(In,Jn,prob_inds);
%     whos
%% do the mapping    
% interim = zeros(tri.size);
tri.size = [tri.size(2) tri.size(1)];
if cx==1
    mapped_img=zeros(tri.size);
else
    mapped_img=zeros([tri.size,dd]);  
    mapped_imgR=zeros([tri.size]);  
    mapped_imgG=zeros([tri.size]); 
    mapped_imgB=zeros([tri.size]); 
    prob_img_cleanR=prob_img_clean(:,:,1);
    prob_img_cleanG=prob_img_clean(:,:,2);
    prob_img_cleanB=prob_img_clean(:,:,3);
end
% whos 
% pause

points1 = probe_mesh_pts;
cont=1;
% for all triangles
if verbose, figure; end

for i=1:tri.count
%     try
        %maketform lahk nadomestim tudi z cp2tform - naredi primerjavo 
        source_triangle_vertices = [probe_pts(tri.indices(i,1),:);probe_pts(tri.indices(i,2),:);probe_pts(tri.indices(i,3),:)];
%         tform = maketform('affine',...  % Make the transformation structure
%                       source_triangle_vertices,...
%                       tri.triangles(i).target);
                  
%          tform = cp2tform(source_triangle_vertices,tri.triangles(i).target,'affine');
%         try
            tform = mycp2tform(source_triangle_vertices,tri.triangles(i).target,'affine');
%         catch
%             non_colin = round((sign(.5-rand(3,2)).*((2-1.5).*rand(3,2)+1.5))+tri.triangles(i).target);
%             tform = mycp2tform(source_triangle_vertices,non_colin,'affine');
%         end

        % poišèimo indekse, ki jih je treba vun vzet iz originala
%         tri.triangles(i)
%         tform.tdata
%         Z = tform.tdata.Tinv'*(tri.triangles(i).coor)';
        Z = tform.'*(tri.triangles(i).coor)';
        probe_x = round(Z(1,:));
        probe_y = round(Z(2,:));
        if cx==1
            mapped_img(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_clean(sub2ind([mn nn],(probe_y),(probe_x)));
            
            if verbose
                imshow(mapped_img,[]);
                hold on
                drawnow 
            end
        else
%             target_size = squeeze([tri.size,dd]);
%             input_size = squeeze([mn nn dd]);
%             target_first = repmat(tri.triangles(i).ref_first,dd,1);
%             target_second = repmat(tri.triangles(i).ref_second,dd,1);
%             traget_third = [ones(size(tri.triangles(i).ref_first));2*ones(size(tri.triangles(i).ref_first));3*ones(size(tri.triangles(i).ref_first))];
%             intput_first = repmat(round(probe_y),dd,1);
%             intput_second = repmat(round(probe_x),dd,1);
%             intput_third = [ones(size(probe_y));2*ones(size(probe_y));3*ones(size(probe_y))];
%             whos 
%             pause
%             mapped_img1(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_clean(sub2ind([mn nn],(probe_y),(probe_x)));
            mapped_imgR(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_cleanR(sub2ind([mn nn],(probe_y),(probe_x)));
            mapped_imgG(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_cleanG(sub2ind([mn nn],(probe_y),(probe_x)));
            mapped_imgB(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_cleanB(sub2ind([mn nn],(probe_y),(probe_x)));
            
            if verbose
                imshow(mapped_imgR,[]);
                hold on
                drawnow 
            end
        end
            
%     catch
% 
%     end
    
    
    
   
end

if cx~=1
    mapped_img(:,:,1)=mapped_imgR;
    mapped_img(:,:,2)=mapped_imgG;
    mapped_img(:,:,3)=mapped_imgB;
%     if verbose
%         figure
%         imshow(uint8(mapped_img),[]);
%     end
end
%% Pomozne funkcije



























