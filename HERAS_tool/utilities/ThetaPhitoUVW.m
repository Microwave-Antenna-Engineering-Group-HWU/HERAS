function [u,v,w] = ThetaPhitoUVW(theta,phi)
% This function computes the angular coordinates theta and phi from the u
% and v coordinates. u=sin(theta)cos(phi), v=sin(theta)sin(phi) and
% w=cos(theta)
u=sin(theta).*cos(phi);
v=sin(theta).*sin(phi);
w=cos(theta);
end