function T_norm = normalise_points(X)
%
% normalise_points determines the homogeneous
% transformation matrix T_norm such that
%
%   X_norm = T_norm*X
%
% defines an X_norm with a mean of 0 and rms
% (root mean squared) distance from the origin of sqrt(2)
 
d = [];
XS = X';
 
s = sum(XS);
x = s(1)/size(XS,1);
y = s(2)/size(XS,1);
z = s(3)/size(XS,1);
 
XSN = XS - [x*ones(size(XS,1),1) y*ones(size(XS,1),1) z*zeros(size(XS,1),1)];
 
for i=1:size(XSN,1),
%    d = [d;sqrt(XSN(i,1)^2+XSN(i,2)^2+XSN(i,3)^2)];
%    d = [d;(XSN(i,1)^2+XSN(i,2)^2+XSN(i,3)^2)];
        d = [d;(XSN(i,1)^2+XSN(i,2)^2)];
end;
%Dm = mean(d);
Dm = sqrt(mean(d));
sf=sqrt(2)/Dm;
for i=1:size(XSN,1),
    XSN(i,1) = XSN(i,1)*sf;
    XSN(i,2) = XSN(i,2)*sf;
    %XSN(i,3) = XSN(i,3)*sf;
end;
 
 
T_norm = [sf 0 -sf*x; 0 sf -sf*y; 0 0 1];