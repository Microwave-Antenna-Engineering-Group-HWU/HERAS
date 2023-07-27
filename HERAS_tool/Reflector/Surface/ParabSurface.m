classdef ParabSurface<POReflectorSurf
    %   Author : Louis Dufour, Salvador Mercader Pellicer                                             
	%   Revisions : 0.1.0 - 25/10/2018  

    properties
        f;      % focal distance
        x_c;    % Offset
        h;      % Clearance
    end

    methods

        function obj = ParabSurface(f, x_c, h)
            % Constructor
            obj.f  = f;
            obj.x_c = x_c;
            obj.h   = h;
        end
        
        function z = getSurf(obj, x, y)
            % This method returns the z coordinate associated to the (x,y)
            % coordinates. All the coordinates are referred to the CS
            % centered in (x_c,0,z_0).
            z_0=obj.x_c^2/(4*obj.f);
            z =(x+obj.x_c).^2/(4*obj.f)+(y).^2/(4*obj.f)-z_0;
        end
        
        function [n_x, n_y, n_z, N] = getNorm(obj, x, y) 
            % This method returns the normal at (x,y) and 
            % N=sqrt(1+dx^2+dy^2). The latter quantity is used to solve
            % surface integrals (Check "Surface integral" page on wikipedia
            % for more information).

            % Compute partial derivative along x and y
            diffzx=(x+obj.x_c)/(2*obj.f);
            diffzy=y/(2*obj.f);

            % Compute N
            N=sqrt(diffzx.^2+diffzy.^2+1);

            % Compute normal components
            n_x=-diffzx./N;
            n_y=-diffzy./N;
            n_z=1./N;
        end
        
        function lines = getLines(~)
            lines = [];
        end
    end
end
