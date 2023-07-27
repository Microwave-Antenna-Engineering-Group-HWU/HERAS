function [k_x,k_y,k_z]=CrossProduct(E_x,E_y,E_z,H_x,H_y,H_z)
% Author : Salvador Mercader Pellicer
k_x=E_y.*H_z-E_z.*H_y;
k_y=E_z.*H_x-E_x.*H_z;
k_z=E_x.*H_y-E_y.*H_x;