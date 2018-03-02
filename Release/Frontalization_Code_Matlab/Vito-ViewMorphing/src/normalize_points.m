% normalize_points - function computing normalization (translation and
% scaling) of the coordinates of the matched points

% xall1 - (n x 1 matrix) x coordinates of points
% yall1 - (n x 1 matrix) y coordinates of points
% xynormalized - (n x 3 matrix) normalized homogeneous coordinates of points 
% T_norm - (3 x 3 matrix) transformation matrix 

function [xynormalized, T_norm] = normalize_points(xall1, yall1)

[how_many_points, x, d] = size(xall1);

mean_x = mean(xall1);
mean_y = mean(yall1);

shifted_xall1 = xall1 - mean_x;
shifted_yall1 = yall1 - mean_y;

sf = sqrt(2)/mean(sqrt(shifted_xall1.^2 + shifted_yall1.^2));

T_norm = [sf 0 -sf*mean_x; 0 sf -sf*mean_y; 0 0 1];

o = ones(how_many_points,1);
xyall1 = [xall1 yall1 o]';

xynormalized = T_norm*xyall1;

xynormalized = xynormalized';