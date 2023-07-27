function [E_x_s,E_y_s,E_z_s]=vectorSph2Cart(V_rho,V_theta,V_phi,theta,phi)
% Authors : Salvador Mercader Pellicer
% Revision - 10/05/23: Francesco Lisi reduced memory requirement
E_x_s=V_rho.*sin(theta).*cos(phi)+V_theta.*cos(theta).*cos(phi)-V_phi.*sin(phi);
E_y_s=V_rho.*sin(theta).*sin(phi)+V_theta.*cos(theta).*sin(phi)+V_phi.*cos(phi);
E_z_s=V_rho.*cos(theta)-V_theta.*sin(theta);
