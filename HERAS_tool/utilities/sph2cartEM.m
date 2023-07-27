function [X,Y,Z] = sph2cartEM(R,TH,PHI)
% This function computes the cartesian coordinates from the cartesian 
% coordinates used in electromagnetic theory 

TH=pi/2-TH;
[X,Y,Z] = sph2cart(PHI,TH,R);

end
