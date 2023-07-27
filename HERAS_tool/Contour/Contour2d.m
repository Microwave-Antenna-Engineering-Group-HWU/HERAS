classdef (Abstract) Contour2d < Contour
% Author : Francesco Lisi
% Revisions : v0.0.1 - 01/03/2023
% -------------------------------------------------------------------------
% This class defines a closed curve in a plane

    properties
        rotation_matrix = eye(2,2);
        N_dimensions = 2;
    end

    methods 

        function Rotate(obj,rotation_angle)
            % This method computes the rotation matrix associated to
            % rotation_angle
            obj.rotation_matrix=[cos(rotation_angle),-sin(rotation_angle);sin(rotation_angle),cos(rotation_angle)];
        end

    end
end