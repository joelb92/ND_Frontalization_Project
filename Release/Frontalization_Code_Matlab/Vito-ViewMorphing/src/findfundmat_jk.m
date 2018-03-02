% findfund:
% calculates a fundamental matrix for the two specified images
%
% im1, im2 - names of two images
% F - (3 x 3 matrix) fundamental matrix
% 
% up' * F * u = 0
% up = [up vp]
% u  = [u v]

function F = findfundmat_jk(lm1,lm2)

u = lm1(:,1);
v = lm1(:,2);
up = lm2(:,1);
vp = lm2(:,2);

% m - no of points
[m,~] = size(u);

o = ones(m,1);

% normalize points
[uv, t1] = normalize_points(u,v);
[upvp, t2] = normalize_points(up,vp);

u = upvp(:,1);
v = upvp(:,2);
up = uv(:,1);
vp = uv(:,2);

% compute A - the equation matrix
A = [u.*up u.*vp u v.*up v.*vp v up vp o];

% the singular value decomposition SVD
[~,~,V] = svd(A,0);

% extract column of the smallest singular value - the last column
smallest = V(:,9);

% reshape the column matrix to 3 x 3 fundamental matrix
F = reshape(smallest,3,3)';

% enforce singularity contraint - 
% fundamental matrix must be singular and of rank 2

[U,D,V] = svd(F,0);
r = D(1,1);
s = D(2,2);

F = U*diag([r s 0])*V';

F = t2'*F*t1;
