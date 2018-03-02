% The function generates an alignment model that can be used with
% the function using it ;).
% 
% PROTOTYPE:
%    amodel = generateCMR_model(mpath)
% 
% Outputs:
%   model ... a structure containing the model that we need for alignemnt;
%             the structure contains the following fields
% 
% Inputs:
%   mpath    ... path to the folder where your model data lives
% 
% 
% Author: Vitomir Struc
% Date: 23 September, 2015
% Copyright: University of Ljubljana, University of Notre Dame, 2015
% 
function amodel = generateCMRmodel(mpath)

% do some path mish-mash
curr_path = cd;
cd(mpath)

%% Load and initialize detection stuff
VJdetector.down_sampling_size = 320; 

VJdetector.fdetector = vision.CascadeObjectDetector();
VJdetector.pdetector = vision.CascadeObjectDetector('haarcascade_profileface.xml');
amodel.VJdetector = VJdetector;

%% Intraface detector

option.face_score = 0.3;
option.min_neighbors = 2;
option.min_face_image_ratio = 0.15;
option.compute_pose = true;
% OpenCV face detector model file
xml_file = 'haarcascade_frontalface_alt2.xml';
% load tracking model
load('TrackingModel-xxsift-v1.10.mat');
% load detection model
load('DetectionModel-xxsift-v1.5.mat');
% create face detector handle
fd_h = vision.CascadeObjectDetector(xml_file);
DM{1}.fd_h = fd_h;
amodel.DM = DM;
amodel.TM = TM;

%% Load and initalize fudicial point stuff
load('CMR_model500X1.mat');
amodel.Ldetector = Ldetector;

%% Generate triangle structure for target shape - dealuney
shape_points = amodel.Ldetector.init_model;
amodel.tri = generate_mapping_data(double(shape_points),0);

%% Load average points
avg_pts = load('lm_avg.mat'); 
amodel.lm_avg = avg_pts.lm_avg;


%% Change path back to initial one
cd(curr_path);

end

%% Function that generates Delaunay triangulation structure needed 
% The function precomputes some data used for aligning shapes based on delauney triangulation.
% It generates a structure that can be used by
% "do_texture_mapping_with_delaunay" to align an annotated probe shape with the
% reference shape used as input for this function.
% 
% Protoype:
%     tri = generate_mapping_data(shape_points)
%     
% Inputs:
%     shape_points ... a Nx2 matrix with x and y coordinates of the target shape
%         
% Outputs
%     tri ... a structure with precomputed data needed for the shape mapping
%      
%     tri.indices ... a Nx3 matrix, where each row defines a triangle  
%                     by indexing entries in shape_points
%     tri.maxs ... a 1x2 matrix with the maximum in the x and y direction of
%                  the shape
%     tri.mins ... a 1x2 matrix with the minimum in the x and y direction of
%                  the shape
%     tri.size ... a 1x2 matrix with the size in the x and y direction of
%                  the shape
%     tri.count ... numer of triangles
%     tri.triangles(i)  ... for the i-th triangle contains the follwoing
%                           fiels
%                  .target ... triangle vertices, 3x2 matrix
%                  .xmin ... minimum in x direction in normalized
%                            coordinates
%                  .xmax ... maximum in x direction in normalized
%                            coordinates
%                  .ymin ... minimum in y direction in normalized
%                            coordinates
%                  .ymax ... maximum in y direction in normalized
%                            coordinates
%                  .ind  ... linear indices of all points in the target
%                            shape belonging to the i-th triangle
%                  .tri_mask ... binary mask of the i-th triangle 
%     
%         
% Author: Vitomir Struc
% Date: 3.9.2014
% Copyright: Fe UL, Luks, 2014

function tri = generate_mapping_data(shape_points,verbose)

    %% init
    if nargin < 2
        verbose = 0;
    end

    %% dummy
    tri=[];
    shape_points = round(shape_points);

    %% generate meta-data
    tri.maxs = [max(shape_points(:,1)), max(shape_points(:,2))]; %maximums
    tri.mins = [min(shape_points(:,1)), min(shape_points(:,2))]; %maximums
    tri.size = [abs(tri.maxs(1)-tri.mins(1))+1,abs(tri.maxs(2)-tri.mins(2))+1]; %size

    %% do triangulation 
    
    %normalize coordinates 
    tri.points = [shape_points(:,1)-tri.mins(1)+1,shape_points(:,2)-tri.mins(2)+1]; 
    
    % do delauney
    dres = DelaunayTri(tri.points);
    tri.indices = dres.Triangulation;%delaunay(shape_points(:,1),shape_points(:,2));
    tri.dres = dres;
    
    %get triangle count
    tri.count = size(tri.indices,1);

    % show points
    if verbose
        figure
        h1 = axes;
        triplot(tri.indices,tri.points(:,1),tri.points(:,2))
        axis([0 tri.size(1) 0 tri.size(2)])
        set(h1, 'Ydir', 'reverse')
        disp('Displaying reference shape. May be shown upside-down.')
        disp('Press button to continue ...')
        pause
    end
    
    %create coordinates of target image
    [In,Jn] = ind2sub([tri.size(2) tri.size(1)],1:(tri.size(1)*tri.size(2)));
    
    
    %get triangle indices
    tri.inds = pointLocation(dres,[Jn;In]');
    tri.inds = sparse(In,Jn,tri.inds);
%     whos 
%     pause
    
%     imshow(reshape(inds,tri.size(2), tri.size(1)),[]);

    %construct triangle data
    for i=1:tri.count
       tri.triangles(i) = get_tri_data(tri,i,verbose);
       if verbose
            disp(sprintf('Fisnished processing triangle %i/%i',i,tri.count)); 
       end
    end



end %function end



%% generate triangle data
function triang = get_tri_data(tri,ind,verbose)
    triang = [];
    
    % to so zej koordinate target trikotnika
    triang.target = [tri.points(tri.indices(ind,1),:);tri.points(tri.indices(ind,2),:);tri.points(tri.indices(ind,3),:)];
    
    % pogledam maximume in minimume target trikotnika
    triang.xmin = min(triang.target(:,1));
    triang.xmax = max(triang.target(:,1));
    triang.ymin = min(triang.target(:,2));
    triang.ymax = max(triang.target(:,2));
    
    [triang.ref_first,triang.ref_second]=find(tri.inds==ind); 

    triang.coor = [triang.ref_second, triang.ref_first,ones(size(triang.ref_first))];
end
































