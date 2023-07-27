function [V_phi,V_theta,V_rho]=vectorCart2Sph(E_x_s,E_y_s,E_z_s,theta,phi)
% Authors : Salvador Mercader Pellicer
V_rho=E_x_s.*sin(theta).*cos(phi)+E_y_s.*sin(theta).*sin(phi)+E_z_s.*cos(theta);
V_theta=E_x_s.*cos(theta).*cos(phi)+E_y_s.*cos(theta).*sin(phi)-E_z_s.*sin(theta);
V_phi=-E_x_s.*sin(phi)+E_y_s.*cos(phi);