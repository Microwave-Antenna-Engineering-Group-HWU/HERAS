function [V_x,V_y,V_z]=sph2cartEMvec(V_r,V_theta,V_phi,theta,phi)
% Authors : Salvador Mercader Pellicer
% Revision - 10/05/23: Francesco Lisi reduced memory requirement
V_x= V_r.*sin(theta).*cos(phi)  +V_theta.*cos(theta).*cos(phi) -V_phi.*sin(phi);
V_y= V_r.*sin(theta).*sin(phi)  +V_theta.*cos(theta).*sin(phi) +V_phi.*cos(phi);
V_z= V_r.*cos(theta)            -V_theta.*sin(theta);