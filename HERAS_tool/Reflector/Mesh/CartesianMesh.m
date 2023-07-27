classdef CartesianMesh<POReflectorMesh
    % CartesianMesh implements a regular cartesian mesh on a reflector
    % with a disk aperture. Can easily be extended to non-disk shape: the
    % mesh is a square, and the points outside of the circular rim are
    % assigned a null surface current
    % Integration is trapz, using Matlab built in function and the parallel
    % processing toolbox
    % Authors : Salvador Mercader Pellicer, Louis Dufour
    % v0.1.0 : 25/10/2018
    
    properties
        x, y;
        Nx, Ny, lambda_0;
        rim;
        Jacobian = 1;
        weights;
    end
    
    methods
        
        function obj = CartesianMesh(Nx, Ny, lambda_0) 
            % Constructor
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.lambda_0 = lambda_0;
            obj.intType = 'cartesian';
        end
        
        function meshNow(obj, rim, lines)
            % This method creates the (x,y) mesh of the reflector
            switch rim.type
                case 'elliptic'
                    % Compute (x,y) mesh
                    obj.rim = rim;
                    x0 = obj.rim.x0;
                    y0 = obj.rim.y0;
                    a1 = obj.rim.a1;
                    a2 = obj.rim.a2;
                    obj.x = linspace(x0-a1,x0+a1,obj.Nx);
                    obj.y = linspace(y0-a2,y0+a2,obj.Ny);
                    [obj.xE,obj.yE] = meshgrid(obj.x,obj.y);
                    obj.xI = obj.xE; obj.yI = obj.yE;
                    
                    % Compute weight matrix for integration
                    w_x_single=(obj.x(end)-obj.x(1))/(2*(obj.Nx-1))*[1, 2*ones(1,obj.Nx-2),1];
                    w_y_single=(obj.y(end)-obj.y(1))/(2*(obj.Ny-1))*[1, 2*ones(1,obj.Ny-2),1];
                    [Wx,Wy]=meshgrid(w_x_single,w_y_single);
                    W=Wx.*Wy;
                    obj.weights=W;

                    % Set isMeshed flag to true
                    obj.isMeshed = true;
                otherwise
                    error('CartesianMesh only works with elliptic rims');
            end
        end
        
        function [x1,x2,J] = getGridPoints(obj)
            x1 = obj.xI;
            x2 = obj.yI;
            J = obj.Jacobian;
        end
       
    end
end