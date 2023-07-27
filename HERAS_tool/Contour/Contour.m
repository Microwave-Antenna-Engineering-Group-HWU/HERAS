classdef Contour < matlab.mixin.Copyable
% Author : Francesco Lisi
% Revisions : v0.0.1 - 01/03/2023
% -------------------------------------------------------------------------
% This abstract class defines a genaral contour in R^(N_dimensions)

    properties (Abstract)
        rotation_matrix;
        N_dimensions;
    end

    properties
        center;             % Center point of the contour
        furthest_point;     % Furthest contour point from center
    end

    methods (Abstract)

        inside_mask = IsInside(obj,points)  % Abstract method that returns the mask of points inside the contour
        out=PlotContour(obj)                % Abstract class to plot the contour

    end

end