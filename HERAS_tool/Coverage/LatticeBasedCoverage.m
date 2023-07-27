classdef LatticeBasedCoverage < Coverage
% Author : Francesco Lisi
% Revisions : 0.1.0 - 06/03/2023
% -------------------------------------------------------------------------
% This class defines the coverage from a contour and a lattice.

    properties
        contour;
        lattice;
    end

    methods
        function obj = LatticeBasedCoverage(x,y,z,ex,ey,ez)
        % Constructor
            
            % Call parent class constructor
            obj@Coverage(x,y,z,ex,ey,ez)
        end

        function SetContour(obj,contour)
        % This method sets the contour of the array
            obj.contour = contour;
        end

        function SetLattice(obj,lattice)
        % This method sets the lattice of the array
            obj.lattice = lattice;
        end

        function CreateCoverage(obj)
            % First compute each u-v point in the coverage CS
            [uv_position,index] = obj.lattice.GeneratePointsInsideContour(obj.contour);

            % Define coverage from u-v points
            obj.DefineCoverageFromUV(uv_position(1,:),uv_position(2,:));

            % Set index
            obj.index = index;

        end

    end
end