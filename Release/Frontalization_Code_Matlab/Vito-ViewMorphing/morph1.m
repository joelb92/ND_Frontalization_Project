function [morphed_im,imwarp_pts] = morph1(im1, im2, im1_pts, im2_pts, tri, warp_frac, dissolve_frac)
% MORPH provides the morphing result of the two input images, which
% represent one kind of intermediate state

% Written by Qiong Wang at University of Pennsylvania
% Oct. 9th, 2013

tri_con = tri.ConnectivityList;

% Check whether the number of control points in both images the same
if size(im1_pts, 1) ~= size(im2_pts)
    error('Number of control points in the two images should be equal!');
end

% Pad images at the largest extent if they are not padded
[nr1, nc1, ~] = size(im1);
[nr2, nc2, ~] = size(im2);
if nr1 == nr2 && nc1 == nc2
    nr = nr1;
    nc = nc1;
    im1_pad = im1;
    im2_pad = im2;
else
    nr = max(nr1, nr2); 
    nc = max(nc1, nc2);
    im1_pad = padarray(im1, [nr-nr1, nc-nc1], 'replicate', 'post');
    im2_pad = padarray(im2, [nr-nr2, nc-nc2], 'replicate', 'post');
end

% Initialize
tri_num    = tri.size(1);
im1_warp   = zeros(nr, nc, 3);
im2_warp   = zeros(nr, nc, 3);
sub_array  = [repmat(1:nc, 1, nr); reshape(repmat(1:nr, nc, 1), [1 nr*nc]); ones(1, nr*nc)];
ps_1       = zeros(3, nr * nc);
ps_2       = zeros(3, nr * nc);

% Intermediate points
imwarp_pts = (1 - warp_frac) * im1_pts + warp_frac * im2_pts;

% Loop for each triangle to calculate the transform matrix
row_ind = pointLocation(tri, sub_array(1:2,:)');

for i = 1 : tri_num
    tf_warp               = [imwarp_pts(tri_con(i, :),:)'; ones(1,3)];
    tf_1                  = [im1_pts(tri_con(i, :),:)'; ones(1,3)] / tf_warp;
    tf_2                  = [im2_pts(tri_con(i, :),:)'; ones(1,3)] / tf_warp;
    ps_1(:, row_ind == i) = pixel_limit(tf_1 * sub_array(:, row_ind == i), nr, nc);
    ps_2(:, row_ind == i) = pixel_limit(tf_2 * sub_array(:, row_ind == i), nr, nc);
end

nan_idx = ~isnan(row_ind);
nan_len = sum(nan_idx);
im1_warp(sub2ind(size(im1_warp), repmat(sub_array(2,~isnan(row_ind)),1,3), repmat(sub_array(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)] )) =...
  im1_pad(sub2ind(size(im1_pad), repmat(ps_1(2,~isnan(row_ind)),1,3), repmat(ps_1(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)]));
im2_warp(sub2ind(size(im2_warp), repmat(sub_array(2,~isnan(row_ind)),1,3), repmat(sub_array(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)] )) =...
  im2_pad(sub2ind(size(im2_pad), repmat(ps_2(2,~isnan(row_ind)),1,3), repmat(ps_2(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)]));

% Dissolve warpped images
morphed_im = (1-dissolve_frac) * im1_warp + dissolve_frac * im2_warp;
morphed_im = uint8(morphed_im);

