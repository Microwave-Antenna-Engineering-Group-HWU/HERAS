classdef CSlocal < POobject
    % CSlocal defines a coordinate system, which is usefull to do
    % coordinate system transformation
    % Reflector and feeds are CSlocal with other properties; this class is
    % to be used when a coordinate system of interest does not correspond
    % to any physical object
    properties
    end
    
    methods
        function obj = CSlocal(x,y,z,ex,ey,ez)
            obj.position = [x;y;z];
            obj.direction = [ex(:)/norm(ex(:)) ey(:)/norm(ey(:)) ez(:)/norm(ez(:))];
        end
        
        function disp2fig(obj, ax, CSref)
            if nargin == 1
                ax = [];
            elseif nargin == 2
                CSref = CSlocal(0,0,0,[1,0,0],[0,1,0],[0,0,1]);
            end
            if isempty(ax)
                ax = figure;
            else
                if isa(ax,'matlab.graphics.axis.Axes')
                    axes(ax);  hold on;
                else
                    figure(ax); hold on;
                end
            end
            obj.plotTriade(ax, CSref);
        end
    end
    
end

