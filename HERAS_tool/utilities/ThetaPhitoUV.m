function [u,v] = ThetaPhitoUV(theta,phi)
% This function computes the angular coordinates theta and phi from the u
% and v coordinates. u=sin(theta)cos(phi) and v=sin(theta)sin(phi)
u=sin(theta).*cos(phi);
v=sin(theta).*sin(phi);

end