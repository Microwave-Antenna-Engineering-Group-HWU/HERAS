function [theta,phi] = UVtoThetaPhi(u,v)
% This function computes the angular coordinates theta and phi from the u,
% v and w coordinates. u=sin(theta)cos(phi), v=sin(theta)sin(phi) and
% w=cos(theta)
z=u+1i*v;
phi=wrapTo2Pi(angle(z));
theta=asin(abs(z));


end