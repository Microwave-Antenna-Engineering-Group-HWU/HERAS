function [theta,phi] = UVWtoThetaPhi(u,v,w)
% This function computes the angular coordinates theta and phi from the u,
% v and w coordinates. u=sin(theta)cos(phi), v=sin(theta)sin(phi) and
% w=cos(theta)

[~,theta,phi]=cart2sphEM(u,v,w);


end