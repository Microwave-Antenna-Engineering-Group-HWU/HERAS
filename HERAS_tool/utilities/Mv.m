function M = Mv(v, theta)
% Rotates along vector v (need not be normalized) by an angle theta (in
% radians)
% Author: Louis Dufour
v = v/norm(v);
x = v(1); y = v(2); z = v(3);
u = [x;y;z];
M = [0 -z y;z 0 -x;-y x 0]*sin(theta)+(eye(3)-u*u')*cos(theta)+u*u';
end
