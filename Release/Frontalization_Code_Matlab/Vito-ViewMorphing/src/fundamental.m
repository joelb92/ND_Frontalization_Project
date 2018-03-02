function F1_norm = fundamental(xx1,xx2)
% 
X1 = [xx1'; ones(size(xx1(:,1)'))];
X2 = [xx2'; ones(size(xx2(:,1)'))];
% 
  
T1_norm = normalise_points(X1);
X1_norm = T1_norm*X1;
T2_norm = normalise_points(X2);
X2_norm = T2_norm*X2;
 
 
u1 = X1_norm(1,:)';
v1 = X1_norm(2,:)';
u2 = X2_norm(1,:)';
v2 = X2_norm(2,:)';
 
A = [u2.*u1 u2.*v1 u2 v2.*u1 v2.*v1 v2 u1 v1 ones(size(u1))];
 
[U S V] = svd(A);
 
Fvec = V(:,end);
F_norm =[Fvec(1,1) Fvec(2,1) Fvec(3,1); ...
   Fvec(4,1) Fvec(5,1) Fvec(6,1); ...
   Fvec(7,1) Fvec(8,1) Fvec(9,1) ];
 
% enforce singularity
[US SS VS] = svd(F_norm);
% 
SS(3,3) = 0;
% 
F1_norm = US*SS*VS';
F1_norm = (T2_norm')*(F1_norm)*(T1_norm);