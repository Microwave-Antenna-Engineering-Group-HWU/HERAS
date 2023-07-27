function [coord_x2,coord_y2,coord_z2] = pos2CSfcn(CS1,CS2,coord_x1,coord_y1,coord_z1)
% This function converts the position vector defined by
% (coord_x1,coord_y1,coord_z1) from the coordinate system expressed by the
% structure CS1 to the one expessed by CS2

coord_size=size(coord_x1);

% Extract coordinate system origin position in a global coordinate system
CS1_origin_position=CS1.position(:);
CS2_origin_position=CS2.position(:);

% Extract coordinate system orientation in a global coordinate system
CS1_orientation=CS1.direction;
CS2_orientation=CS2.direction;

% Create position array in CS1
points_1=[coord_x1(:).';coord_y1(:).';coord_z1(:).'];

% Conversion of the points from CS1 to CS2
points_2 = CS2_orientation'*(CS1_orientation*points_1+CS1_origin_position-CS2_origin_position);

% Output extraction
coord_x2 = reshape(points_2(1,:),coord_size);
coord_y2 = reshape(points_2(2,:),coord_size);
coord_z2 = reshape(points_2(3,:),coord_size);

end