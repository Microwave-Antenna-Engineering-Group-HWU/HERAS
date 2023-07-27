classdef (Abstract) POReflectorMesh < matlab.mixin.Copyable
    % POReflectorMesh defines the mesh and integration methods over a
    % reflector. This abstract class must be implemented by any mesh
    % subclass.
    %
    % Properties:
    %   xE, yE, zE  - the coordinates (in the reflector coordinate system)
    %                 of the mesh vertices 
    %   xI, yI, zI  - the coordinates (in the reflector coordinate system)
    %                 of the mesh integratio points. The mesh integration 
    %                 points may (trapz) or may not (Gauss-Legendre) be the
    %                 same
    %   n_x, n_y, n_z, N  - the normals and Jacobian(?) at the integration
    %                points
    %   isMeshed   - false at object creation, take notes if the mesh has 
    %                been meshed or not yet 
    %   displayOptions - structure with display options, to be defined in
    %                the subclass
    properties
        % Points defining the elements
        xE, yE,zE
        % Points where values must be evaluated
        xI,yI, zI;
        % Normals at the integration points xI, yI
        n_x, n_y, n_z, N;
        % 
        isMeshed = false;
        displayOptions;
        intType;
    end
    
    methods (Abstract)
        meshNow(obj, rim, lines);
        %integrateFF(obj);
        %displayMesh(obj);
    end
    
    methods
        function [x,y,z] = getPointsE(obj)
            x = obj.xE;
            y = obj.yE;
            z = obj.zE;
        end
        
        function [x,y,z] = getPointsI(obj)
            x = obj.xI;
            y = obj.yI;
            z = obj.zI;
        end
        
        function yn = ismeshed(obj)
            yn = obj.isMeshed;
        end
        
        function setZE(obj,z)
            obj.zE = z;
        end
        
        function setZI(obj,z)
            obj.zI = z;
        end
        
        function setNormals(obj, n_x, n_y, n_z, N)
            obj.n_x = n_x;
            obj.n_y = n_y;
            obj.n_z = n_z;
            obj.N = N;
        end
        
        function [n_x, n_y, n_z, N] = getNormals(obj)
            n_x = obj.n_x;
            n_y = obj.n_y;
            n_z = obj.n_z;
            N = obj.N;
        end
        
    end
end