function [H1,H2] = PreWarp1(F)

[U S V] = svd(F);

 % Find eigen vectors of F
 e1 = V(1:3,3) / V(3,3);
 d1 = [-e1(2); e1(1); 0];
 
 % find eigen vector of F'
 e2 = U(1:3,3) /U(3,3);

%mine
%[vs,ds] = eig(F);
%min_col = find(diag(ds)==min(diag(ds)));
%e1 = vs(:,min_col);
%d1 = [-e1(2); e1(1); 0];

%[vs,ds] = eig(F');
%min_col = find(diag(ds)==min(diag(ds)));
%e2 = vs(:,min_col);
%%%%%%

Fd1 = F * d1;
Fd1 = Fd1 / sqrt(sum(Fd1.^2));

d2 = [-Fd1(2); Fd1(1); 0];
% calculate angle theta
theta1 = atan(e1(3)/(d1(2)*e1(1) - d1(1)*e1(2)));
theta2 = atan(e2(3)/(d2(2)*e2(1) - d2(1)*e2(2)));
% calculate rotation matrices
RDT1 = Rotation1(d1, theta1);
RDT1N = Rotation1(d1, -1*theta1);
RDT2 = Rotation1(d2, theta2);
% calculate rotated epipolar vectors
eRDT1 = RDT1 * e1;
eRDT2 = RDT2 * e2;
% calculate angle phi
phi1 = -atan(eRDT1(2)/eRDT1(1));
phi2 = -atan(eRDT2(2)/eRDT2(1));
% calculate the phi rotation
RZP1 = RotationPhi(phi1);
RZP1N = RotationPhi(-phi1);
RZP2 = RotationPhi(phi2);

FR = RZP2 * RDT2 * F * RDT1N * RZP1N;

FR = FR / FR(3,2);

T = [1 0 0; 0 -FR(2,3) -FR(3,3); 0 0 1];

H1 = RZP1 * RDT1;
H2 = T * RZP2 * RDT2;