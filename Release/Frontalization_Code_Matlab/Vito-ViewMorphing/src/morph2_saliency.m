function [morphed_im,morphed_pts] = morph2_saliency(im1, im2, lm1, lm2, varargin)
% MORPH provides the morphing result of the two input images, which
% represent one kind of intermediate state

% Required inputs
p = inputParser;
addRequired(p,'im1',@isnumeric)
addRequired(p,'im2',@isnumeric)
addRequired(p,'lm1',@isnumeric)
addRequired(p,'lm2',@isnumeric)

% Optional name-value pairs
addParameter(p,'warp_frac',0.5,@(x)x>=0&&x<=1)
addParameter(p,'dissolve_frac',0.5,@(x)x>=0&&x<=1)
addParameter(p,'yaw_angle',[],@isnumeric)

% Parse inputs
parse(p, im1, im2, lm1, lm2, varargin{:});
im1 = p.Results.im1;
im2 = p.Results.im2;
lm1 = p.Results.lm1;
lm2 = p.Results.lm2;
warp_frac = p.Results.warp_frac;
dissolve_frac = p.Results.dissolve_frac;
yaw_angle = p.Results.yaw_angle;


% Check whether the number of control points in both images the same
if size(lm1, 1) ~= size(lm2)
    error('Number of control points in the two images should be equal!');
end

morphed_pts = (1 - warp_frac) * lm1 + warp_frac * lm2;

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
    im1_pad = padarray(im1, [nr-nr1, nc-nc1], 'replicate', 'post');%imresize(im1, [nr nc]);
%     lm1 = lm1.*repmat([nc/nc1 nr/nr1],length(lm1),1);
    im2_pad = padarray(im2, [nr-nr2, nc-nc2], 'replicate', 'post');%imresize(im2, [nr nc]);
%     lm2 = lm2.*repmat([nc/nc2 nr/nr2],length(lm2),1);
end

points_avg = lm1*(1-warp_frac) + lm2*warp_frac;
[points_avg,idx] = remove_duplicate_points(points_avg);
tri = delaunayTriangulation(points_avg);
tri_con = tri.ConnectivityList;
lm1 = lm1(idx,:);
lm2 = lm2(idx,:);

% Initialize
tri_num    = tri.size(1);
if size(im1,3) == 3
    im1_warp   = zeros(nr, nc, 3);
else
    im1_warp   = zeros(nr, nc);
end
if size(im2,3) == 3
    im2_warp   = zeros(nr, nc, 3);
else
    im2_warp   = zeros(nr, nc);
end
sub_array  = [repmat(1:nc, 1, nr); reshape(repmat(1:nr, nc, 1), [1 nr*nc]); ones(1, nr*nc)];
ps_1       = zeros(3, nr * nc);
ps_2       = zeros(3, nr * nc);

% Intermediate points
imwarp_pts = (1 - warp_frac) * lm1 + warp_frac * lm2;

% Loop for each triangle to calculate the transform matrix
row_ind = pointLocation(tri, sub_array(1:2,:)');

for i = 1 : tri_num
    tf_warp               = [imwarp_pts(tri_con(i, :),:)'; ones(1,3)];
    tf_1                  = [lm1(tri_con(i, :),:)'; ones(1,3)] / tf_warp;
    tf_2                  = [lm2(tri_con(i, :),:)'; ones(1,3)] / tf_warp;
    ps_1(:, row_ind == i) = pixel_limit(tf_1 * sub_array(:, row_ind == i), nr, nc);
    ps_2(:, row_ind == i) = pixel_limit(tf_2 * sub_array(:, row_ind == i), nr, nc);
end

nan_idx = ~isnan(row_ind);
nan_len = sum(nan_idx);
if size(im1_warp,3) == 3
    im1_warp(sub2ind(size(im1_warp), repmat(sub_array(2,~isnan(row_ind)),1,3), repmat(sub_array(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)] )) =...
      im1_pad(sub2ind(size(im1_pad), repmat(ps_1(2,~isnan(row_ind)),1,3), repmat(ps_1(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)]));
    im2_warp(sub2ind(size(im2_warp), repmat(sub_array(2,~isnan(row_ind)),1,3), repmat(sub_array(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)] )) =...
      im2_pad(sub2ind(size(im2_pad), repmat(ps_2(2,~isnan(row_ind)),1,3), repmat(ps_2(1,~isnan(row_ind)),1,3), [ones(1,nan_len) 2*ones(1,nan_len) 3*ones(1,nan_len)]));
else
    im1_warp(sub2ind(size(im1_warp), sub_array(2,~isnan(row_ind)), sub_array(1,~isnan(row_ind)) )) =...
      im1_pad(sub2ind(size(im1_pad), ps_1(2,~isnan(row_ind)), ps_1(1,~isnan(row_ind)) ));
    im2_warp(sub2ind(size(im2_warp), sub_array(2,~isnan(row_ind)), sub_array(1,~isnan(row_ind)) )) =...
      im2_pad(sub2ind(size(im2_pad), ps_2(2,~isnan(row_ind)), ps_2(1,~isnan(row_ind)) ));
end

% Remove all-zero rows and columns
fst_nzero_row = find(any(im1_warp(:,:,1),2),1);
fst_nzero_col = find(any(im1_warp(:,:,1),1),1);
lst_nzero_row = find(any(im1_warp(:,:,1),2)); lst_nzero_row = lst_nzero_row(end);
lst_nzero_col = find(any(im1_warp(:,:,1),1)); lst_nzero_col = lst_nzero_col(end);
im1_warp = im1_warp(fst_nzero_row:lst_nzero_row,fst_nzero_col:lst_nzero_col,:);
im2_warp = im2_warp(fst_nzero_row:lst_nzero_row,fst_nzero_col:lst_nzero_col,:);
imwarp_pts = imwarp_pts - repmat([fst_nzero_col,fst_nzero_row],length(imwarp_pts),1) + 1;
morphed_pts = morphed_pts - repmat([fst_nzero_col,fst_nzero_row],length(morphed_pts),1) + 1;

% Dissolve warpped images
if isempty(yaw_angle) % linear cross dissolve
    
    morphed_im = (1-dissolve_frac) * im1_warp + dissolve_frac * im2_warp;
        
else % saliency preserving and contrast enhancement
   
    mask_eyer = zeros(size(im1_warp,1),size(im1_warp,2));
    mask_eyer(sub2ind(size(mask_eyer),round(imwarp_pts(43:48,2)),round(imwarp_pts(43:48,1)))) = 1;
    mask_eyer = bwconvhull(mask_eyer);
    mask_eyel = zeros(size(im1_warp,1),size(im1_warp,2));
    mask_eyel(sub2ind(size(mask_eyel),round(imwarp_pts(37:42,2)),round(imwarp_pts(37:42,1)))) = 1;
    mask_eyel = bwconvhull(mask_eyel);
    
    mask = repmat(linspace(0,1,size(im1_warp,2)),size(im1_warp,1),1,size(im1_warp,3));
    if yaw_angle < 0
        mask(repmat(mask_eyer,1,1,size(im1_warp,3))) = 1;
        mask(repmat(mask_eyel,1,1,size(im1_warp,3))) = 1;
        mask = imfilter(mask, ones(20,20)/20^2,'replicate');
        
    else
        mask(repmat(mask_eyer,1,1,size(im1_warp,3))) = 0;
        mask(repmat(mask_eyel,1,1,size(im1_warp,3))) = 0;
        mask = imfilter(mask, ones(20,20)/20^2,'replicate');
        mask = (1-mask);
    end
    
    if 0
        
        morphed_im = mask.*im1_warp + (1-mask).*im2_warp;
        
    else

        gama = 10;
        mask1 = mask.^gama./(mask.^gama + (1-mask).^gama); %Eqn. 20 (salience matte for 1st warp)
        mask2 = (1-mask).^gama./(mask.^gama + (1-mask).^gama); %Eqn. 20 (salience matte for 2nd warp)
        morphed_im = mask1.*im1_warp + mask2.*im2_warp;

        % Contrast enhancement (Eqn. 5 - Eqn. 9)
        mu1 = sum(sum(mask.*im1_warp))./sum(sum(mask));
        mu2 = sum(sum((1-mask).*im2_warp))./sum(sum(1-mask));
        sigma1 = sqrt(sum(sum(mask.*(im1_warp-repmat(mu1,size(im1_warp,1),size(im1_warp,2),1)).^2))./sum(sum(mask)));
        sigma2 = sqrt(sum(sum((1-mask).*(im2_warp-repmat(mu2,size(im2_warp,1),size(im2_warp,2),1)).^2))./sum(sum(1-mask)));
        sigma12 = sum(sum(sqrt(mask.*(1-mask)).*(im1_warp-repmat(mu1,size(im1_warp,1),size(im1_warp,2),1)).*(im2_warp-repmat(mu2,size(im1_warp,1),size(im1_warp,2),1))))./sum(sum(sqrt(mask.*(1-mask))));
        sigma = sqrt((mask.*repmat(sigma1,size(im1_warp,1),size(im1_warp,2),1)).^2 + ...
               (1-mask).*repmat(sigma2,size(im1_warp,1),size(im1_warp,2),1).^2 + ...
               2*mask.*(1-mask).*repmat(sigma12,size(im1_warp,1),size(im1_warp,2),1));
        mu = mask.*repmat(mu1,size(im1_warp,1),size(im1_warp,2),1) + ...
         (1-mask).*repmat(mu2,size(im2_warp,1),size(im2_warp,2),1);
        sigma_ = mask.*repmat(sigma1,size(im1_warp,1),size(im1_warp,2),1) + ...
             (1-mask).*repmat(sigma2,size(im2_warp,1),size(im2_warp,2),1);
        tau = 1.05;
        morphed_im = uint8(tau*sigma_./sigma.*(double(morphed_im) - mu) + mu);
        
    end
end

morphed_im = uint8(morphed_im);
