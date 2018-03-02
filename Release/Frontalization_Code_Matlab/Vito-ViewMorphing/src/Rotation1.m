function R=Rotation1(n,theta)

% r' = rcos(theta) + n(n r)(1 - cos(theta) + (r X n) sin(theta)

c = cos(theta);
s = sin(theta);
t = 1 - cos(theta);
x = n(1);
y = n(2);


R = [ x*x+(1-x*x)*c     x*y*t           y*s; ...
      x*y*t             y*y+(1-y*y)*c  -x*s; ...
      -y*s              x*s             c   ];