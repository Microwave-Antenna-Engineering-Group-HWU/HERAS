classdef GeneralEllipseReflector<POReflector
    % This class describes a general elliptic rim reflector, which
    % specified semi major axis and center
    %   Author : Louis Dufour                                             
	%   Revisions : 0.2.0 - 19/11/2018  
    properties
    end
    methods
        function obj = GeneralEllipseReflector(a1,a2,x0,y0,x,y,z)
            %x_c=h+d/2;
            obj.rim.type = 'elliptic';
            obj.rim.a1 = a1;
            obj.rim.a2 = a2;
            obj.rim.x0 = x0;
            obj.rim.y0 = y0;
            obj.position = [x;y;z];
            obj.direction = eye(3);
        end
        
        function inside_mask = InsideRimMask(obj)
            [x,y,~] = obj.getMeshPointsI('local');
            inside_mask=((x-obj.rim.x0).^2/obj.rim.a1^2+(y-obj.rim.y0).^2/obj.rim.a2^2)<=1;
        end

        function out = nullSurfValue(obj,in)
            out = in;
            [x,y,~] = obj.getMeshPointsI('local');
            out(((x-obj.rim.x0).^2/obj.rim.a1^2+(y-obj.rim.y0).^2/obj.rim.a2^2)>1.001)=0;
        end
        function out = nanSurfValue(obj,in)
            out = in;
            [x,y,~] = obj.getMeshPointsE('local');
            out(((x-obj.rim.x0).^2/obj.rim.a1^2+(y-obj.rim.y0).^2/obj.rim.a2^2)>1.001)=nan;
        end
        
        function [x0, y0, z0] = getRimCenter(obj, option)
             x0 = obj.rim.x0;
             y0 = obj.rim.y0;
             z0 = obj.surface.getSurf(x0,y0);
             if strcmp(option, 'global')
                CSglobal  = CSlocal(0,0,0,[1 0 0],[0 1 0],[0 0 1]);
                [x0,y0,z0] = CSglobal.pos2CS(x0,y0,z0,obj);
             end
         end
        
    end
end