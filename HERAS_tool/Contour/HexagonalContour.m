classdef HexagonalContour < Contour2d
% Author : Francesco Lisi
% Revisions : v0.0.1 - 01/03/2023
% -------------------------------------------------------------------------
% This class defines a hexagonal contour in the plane.

    properties
        edge_length;
    end

    methods

        function obj = HexagonalContour(edge_length,varargin)
        % Constructor 
            
            % The second input corresponds to the contour center, which
            % correponds to the origin as default.
            obj.edge_length = edge_length;
            if nargin==2
                obj.center = varargin{1};
            else
                obj.center = zeros(2,1);
            end

            % In an hexagon the furthest point distance corresponds to the
            % edge length
            obj.furthest_point=obj.edge_length;
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
            % edge length and equilateral triangle height
            l=obj.edge_length;
            h=l*sin(pi/3);

            points=obj.rotation_matrix.'*(points-obj.center);
            x=points(1,:);
            y=points(2,:);
            inside_mask=( abs(y) <= h ) & ( abs(y-(2*h/l)*x) <= (2*h) ) & ( abs(y+(2*h/l)*x) <= (2*h) );
            if nargin==3
                inside_mask=reshape(inside_mask,size(varargin{1}));
            end        
        end

        function varargout = PlotContour(obj,varargin)
            phi=0:pi/3:(2*pi);
            points=obj.rotation_matrix*[obj.edge_length*cos(phi);obj.edge_length*sin(phi)]+obj.center;
            varargout{:}=plot(points(1,:),points(2,:),'r',varargin{:},'LineWidth',1.5);
        end

    end
end