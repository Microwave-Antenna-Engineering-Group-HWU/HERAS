classdef RectangularContour < Contour2d
% Author : Francesco Lisi
% Revisions : v0.0.1 - 01/03/2023
% -------------------------------------------------------------------------
% This class defines a rectangular contour in the plane.

    properties
        l_x;
        l_y;
    end

    methods

        function obj = RectangularContour(l_x,varargin)
        % Constructor 
            obj.l_x = l_x;
            switch nargin
                case 1
                    obj.l_y=l_x;
                    obj.center = zeros(2,1);
                case 2
                    if size(varargin{1},1)==1
                        obj.l_y=varargin{1};
                        obj.center = zeros(2,1);
                    else
                        obj.center=varargin{1};
                        obj.l_y=obj.l_x;
                    end
                case 3
                    obj.l_y=varargin{1};
                    obj.center=varargin{2};
            end
            obj.furthest_point=max(obj.l_x,obj.l_y)/2;
        end

        function inside_mask = IsInside(obj,varargin)
        % This function checks if the points in points are inside the
        % rectangular region, and returns a mask, where the value 1 
        % correponds to a point inside the region.
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
            points=obj.rotation_matrix.'*(points-obj.center);
            x=points(1,:);
            y=points(2,:);
            inside_mask=(abs(x)<=(obj.l_x/2))&(abs(y)<=(obj.l_y/2));
            if nargin==3
                inside_mask=reshape(inside_mask,size(varargin{1}));
            end
        end

        function varargout = PlotContour(obj,varargin)
            phi=pi/4:pi/2:(2*pi+pi/4);
            vertex_matrix=sqrt(2)*[(obj.l_x/2)*cos(phi);(obj.l_y/2)*sin(phi)];
            points=obj.rotation_matrix*vertex_matrix+obj.center;
            varargout{:}=plot(points(1,:),points(2,:),'r',varargin{:},'LineWidth',1.5);
        end

    end
end