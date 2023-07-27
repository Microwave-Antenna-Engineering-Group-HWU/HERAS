function [V_r,V_theta,V_phi]=cart2sphEMvec(V_x,V_y,V_z,theta,phi)
% Authors : Salvador Mercader Pellicer
V_r     =+V_x.*sin(theta).*cos(phi) +V_y.*sin(theta).*sin(phi) +V_z.*cos(theta);
V_theta =+V_x.*cos(theta).*cos(phi) +V_y.*cos(theta).*sin(phi) -V_z.*sin(theta);
V_phi   =-V_x.*sin(phi)             +V_y.*cos(phi);