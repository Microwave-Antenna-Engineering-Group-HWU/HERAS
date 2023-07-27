classdef CircleReflector<POReflector
    % This is a parabolic reflector, centered on the far field CS, and
    % alligned with the far field CS. The far field CS is also the global
    % CS
    % The rim is circular, of diameter d. The type of rim is 'disk'
    %   Author : Louis Dufour                                             
	%   Revisions : 0.1.0 - 25/10/2018  
    properties
    end
    methods
        function obj = CircleReflector(d,x,y,z)
            %x_c=h+d/2;
            %obj.d = d;
            %obj.rim = @(phi) d/2;
            obj.rim.type = 'elliptic';
            obj.rim.a1 = d/2;
            obj.rim.a2 = d/2;
            obj.rim.x0 = 0;
            obj.rim.y0 = 0;
            obj.position = [x;y;z];
            obj.direction = eye(3);
        end
        
        function inside_mask = InsideRimMask(obj)
            [x,y,~] = obj.getMeshPointsI('local');
            radius = obj.rim.a1;
            inside_mask=sqrt(x.^2+y.^2)<=radius;
        end

        function out = nullSurfValue(obj,in)
            tol=1e-6;
            out = in;
            [x,y,~] = obj.getMeshPointsI('local');
            d = obj.rim.a1*2;
            out(sqrt(x.^2+y.^2)>(d/2+tol))=0;
        end
        function out = nanSurfValue(obj,in)
            tol=1e-6;
            out = in;
            [x,y,~] = obj.getMeshPointsE('local');
            d = obj.rim.a1*2;
            out(sqrt(x.^2+y.^2)>(d/2+tol))=nan;
        end
        
         function [x0, y0, z0] = getRimCenter(obj, option)
             x0 = 0;
             y0 = 0;
             z0 = obj.surface.getSurf(x0,y0);
             if strcmp(option, 'global')
                CSglobal  = CSlocal(0,0,0,[1 0 0],[0 1 0],[0 0 1]);
                [x0,y0,z0] = CSglobal.pos2CS(x0,y0,z0,obj);
             end
         end
        
    end
end