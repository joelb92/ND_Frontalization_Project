function [ H1, H2 ] = computeH1H2( F )

[V,D] = eig(F);
% i = (diag(D) >= -0.001) & (diag(D) <= 0.001); %original
i = abs(diag(D)) == min(abs(diag(D))); %jk
e0 = V(:,i);

[V,D] = eig(F');
% i = (diag(D) >= -0.001) & (diag(D) <= 0.001); %original
i = abs(diag(D)) == min(abs(diag(D))); %jk
e1 = V(:,i);

d0 = [0;e0(1);0];%[-e0(2);e0(1);0];

Fd0 = F*d0;
d1 = [0;Fd0(1);0];%[-Fd0(2);Fd0(1);0];

theta0 = atan(e0(3)/(d0(2)*e0(1) - d0(1)*e0(2)));
theta1 = atan(e1(3)/(d1(2)*e1(1) - d1(1)*e1(2)));

c = cos(theta0);
s = sin(theta0);
t = 1-cos(theta0);
x = d0(1);
y = d0(2);

R0 = [t*x*x + c   t*x*y       s*y;...
      t*x*y       t*y*y + c   -s*x;...
      -s*y         s*x        c];

c = cos(theta1);
s = sin(theta1);
t = 1-cos(theta1);
x = d1(1);
y = d1(2);

R1 = [t*x*x + c   t*x*y       s*y;...
      t*x*y       t*y*y + c   -s*x;...
      -s*y         s*x        c];


Re0 = R0*e0;
Re1 = R1*e1;

phi0 = -atan(Re0(2)/Re0(1));
phi1 = -atan(Re1(2)/Re1(1));

Rphi0 = [cos(phi0) -sin(phi0) 0;sin(phi0) cos(phi0) 0;0 0 1];
Rphi1 = [cos(phi1) -sin(phi1) 0;sin(phi1) cos(phi1) 0;0 0 1];

H1 = Rphi0*R0;
H2 = Rphi1*R1;

% H1 = R0;
% H2 = R1;

Fw = H2*F*H1';

if ~( all((2*abs(diag(H1))) >= sum(abs(H1),2)) && all((2*abs(diag(H2))) >= sum(abs(H2),2)) )
    H1 = eye(3);
    H2 = H1;
end
