classdef EllipticalContour < Contour2d
% Author : Francesco Lisi
% Revisions : v0.0.1 - 01/03/2023
% -------------------------------------------------------------------------
% This class defines an elliptical contour in the plane.

    properties
        x_half_axis;
        y_half_axis;
    end

    methods

        function obj = EllipticalContour(varargin)
        % Constructor 

            % If nargin==1 then the only input corresponds to the radius of
            % a circular contour centered in the origin.
            % If nargin==2  then the two inputs correspond to the minor and 
            % major axis of an elliptical contour centered in the origin.
            % If nargin==3  then the three inputs correspond to the minor 
            % axis, major axis, and center of an elliptical contour.
            switch nargin
                case 1
                    obj.x_half_axis=varargin{1};
                    obj.y_half_axis=varargin{1};
                    obj.center=zeros(2,1);
                case 2
                    if numel(varargin{2})==1
                        obj.x_half_axis=varargin{1};
                        obj.y_half_axis=varargin{2};
                        obj.center=zeros(2,1);
                    else
                        obj.x_half_axis=varargin{1};
                        obj.y_half_axis=varargin{1};
                        obj.center=varargin{2};
                    end
                case 3
                    obj.x_half_axis=varargin{1};
                    obj.y_half_axis=varargin{2};
                    obj.center=varargin{3};
            end

            % COmpute furthest_point
            obj.furthest_point=max(obj.x_half_axis,obj.y_half_axis);

        end

        function inside_mask = IsInside(obj,varargin)
        % This function checks if the points in points are inside the
        % circular region, and returns a mask, where the value 1 correponds
        % to a point inside the region.

            % If the user passes one input this must be a 2xN vector
            % containing the x and y coordinate of each point. If the user
            % passes two inputs, these must be two equal size matrixes 
            % containing the x and y coordinates of each point.
            switch nargin
                case 2
                    points=varargin{1};
                case 3
                    if isequal(size(varargin{1}),size(varargin{2}))
                        N=numel(varargin{1});
                        points=zeros(2,N);
                        points(1,:)=reshape(varargin{1},[1,N]);
                        points(2,:)=reshape(varargin{2},[1,N]);
                    else
                        error('Wrong input size: the two inputs should be the same size.');
                    end
            end

            % Adjust for the rotation of the contour
            points=obj.rotation_matrix.'*(points-obj.center);

            % Extract x and y coordinates
            x=points(1,:);
            y=points(2,:);

            % COmpute inside_mask
            inside_mask=((x/obj.x_half_axis).^2+(y/obj.y_half_axis).^2)<=1;

            % Reshape inside_mask
            if nargin==3
                inside_mask=reshape(inside_mask,size(varargin{1}));
            end
        end

        function varargout = PlotContour(obj,varargin)
            % Plot elliptical contour
            phi=0:(pi/2^10):(2*pi);
            vertex_matrix=[obj.x_half_axis*cos(phi);obj.y_half_axis*sin(phi)];
            points=obj.rotation_matrix*vertex_matrix+obj.center;
            varargout{:}=plot(points(1,:),points(2,:),'r',varargin{:},'LineWidth',1.5);
        end

    end
end