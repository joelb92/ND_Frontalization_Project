%% The function performs texture mapping based on delaunay triangualtion. 
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

function mapped_img = do_3Dtexture_mapping_with_delaunay(probe_mesh_pts, prob_image, img,tri,P,verbose)

%% init
if nargin<4
    verbose = 0;
end

%% clean probe image - remove all parts of image except for the face covered by the landmarks
max_x = max(probe_mesh_pts(:,1));
min_x = min(probe_mesh_pts(:,1));
max_y = max(probe_mesh_pts(:,2));
min_y = min(probe_mesh_pts(:,2));
size_y = max_y-min_y+1;
size_x = max_x-min_x+1;

[ax,bx,cx]=size(prob_image);
prob_img_clean=prob_image;
if cx==1
%     prob_img_clean = prob_image(min_y:max_y,min_x:max_x);
    [mn nn] = size(prob_img_clean);
else
%     prob_img_clean = prob_image(min_y:max_y,min_x:max_x,:);
    [mn nn dd] = size(prob_img_clean);
end

% normalize probe shape coordinates and find triangle points
probe_mesh_pts=round(probe_mesh_pts);
% probe_pts = [probe_mesh_pts(:,1)-min_x+1,probe_mesh_pts(:,2)-min_y+1];
probe_pts=probe_mesh_pts;
%% Initialize output  
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


%% do the mapping
points1 = probe_mesh_pts;
cont=1;
for i=1:tri.count % for all triangles
        Z = P*(tri.triangles(i).coor)';
        probe_x = round(Z(1,:));
        probe_y = round(Z(2,:));
        
%         for k=1:size(tri.triangles(i).coor,1)
%             img(tri.triangles(i).coor(k,2),tri.triangles(i).coor(k,1))=0;
%         end
%         
%         
%         for k=1:length(probe_x)
%             prob_image(probe_y(k),probe_x(k),:)=0;
%         end
%         figure
%         imshow(prob_image,[])
%         figure(1)
%         imshow(img,[])
%         drawnow
%         pause(0.01)
        
        
%         source_triangle_vertices = [probe_pts(tri.indices(i,1),:);probe_pts(tri.indices(i,2),:);probe_pts(tri.indices(i,3),:)];
        if cx==1
            mapped_img(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_clean(sub2ind([mn nn],(probe_y),(probe_x)));           
            if verbose
                figure(1)
                imshow(mapped_img,[]);
                drawnow 
            end
        else
            mapped_imgR(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_cleanR(sub2ind([mn nn],(probe_y),(probe_x)));
            mapped_imgG(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_cleanG(sub2ind([mn nn],(probe_y),(probe_x)));
            mapped_imgB(sub2ind(tri.size,tri.triangles(i).ref_first,tri.triangles(i).ref_second)) = prob_img_cleanB(sub2ind([mn nn],(probe_y),(probe_x)));
            
            if verbose
                figure(2)
                imshow(mapped_imgR,[]);
                drawnow 
            end
        end  
end
% figure
% imshow(prob_image,[])
% figure
% imshow(img,[])
% pause

if cx~=1
    mapped_img(:,:,1)=mapped_imgR;
    mapped_img(:,:,2)=mapped_imgG;
    mapped_img(:,:,3)=mapped_imgB;
end
end