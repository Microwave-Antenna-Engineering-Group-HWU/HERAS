function [R,TH,PHI] = cart2sphEM(X,Y,Z)
% This function computes the spherical coordinates used in electromagnetic
% theory from the cartesian coordinates

[PHI,TH,R] = cart2sph(X,Y,Z);
TH=pi/2-TH;

end
