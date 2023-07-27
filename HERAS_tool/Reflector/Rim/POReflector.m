classdef (Abstract) POReflector < POobject
    % POReflector describes a general reflector. It defines the interfaces to be 
    %   implemented by any feed to integrate seemlessly with the other parts
    %   of the PO tool, and provide basic methods common to all reflectors.
    %   It serves as a container for the mesh and surface of the reflector,
    %   and the only specific thing it implements is the definition of the
    %   rim.
    %   It inherits from the abstract class POobject
    %   Author : Louis Dufour                                             
    %   Revisions : 0.1.0 - 25/10/2018
    %   Revisions : 0.1.1 - 06/11/2018 - integrationFF moved to dedicated
    %   function file, and 'local' and 'global' introduced for get points
    %
    %   Properties:
    %       surface     - a POReflectorSurf object
    %       mesh        - a POReflectorMesh object
    %       rim         - a structure that says the type of rim, and its parameters
    %
    %  Methods (Abstract)
    %       nullSurfValue  - replace with 0 all values that are outside of
    %                        the reflector rim. Points corresponf to
    %                        integration points.
    %       nanSurfValue   - replace with nan all values that are outside of
    %                        the reflector rim. Points correspond to mesh
    %                        vertices.
    %
    %  Methods
    %       setMesh     - assigns a mesh to the reflector
    %       setSurface  - assigns a surface to the reflector
    %       getMeshPointsI  - get the integration points of the mesh
    %                         associated with the reflector
    %       getMeshPointsE  - get the vertices points of the mesh
    %                         associated with the reflector
    %       getMeshNormals  - get the normals of the reflector surface at
    %                         the mesh integration points
    %       depthReflector - once the mesh and the surface is assigned,
    %                        calculates the z positions of both vertices 
    %                        and mesh integration points, and the normals 
    %                        at the mesh integration points
    %       integrateFF    - integrate the far field given the surface
    %                        currents
    %       displayReflector  - display the reflector, either in the figure
    %                       given in parameters or in a new figure
    %   See POobject
    
    properties
        surface; % a POReflectorSurf object defining the reflector surface geometry
        mesh;    % a POReflectorMesh object defining the mesh and the integration methods
        rim;     % a struct defining the reflector rim. The content of the structure changes with each
                 % reflector type, but must at least contain the field
                 % 'type'
    end
    
    methods (Abstract)
        mask = InsideRimMask(obj);
        out = nullSurfValue(obj,in);
        out = nanSurfValue(obj,in);
        [x0, y0, z0] = getRimCenter(obj); % return the center of the reflector rim in the local or global CS
    end
    
    methods (Access = protected)
      function cp = copyElement(obj)
         % Shallow copy object
         cp = copyElement@matlab.mixin.Copyable(obj);
         % Deep copy for the mesh
         obj.mesh = obj.mesh.copy();
         % Deep copy for the surface
         obj.surface = obj.surface.copy();
      end
   end
    
    methods
        function setMesh(obj, mesh)
            % setMesh assigns a mesh onto the reflector
            %   obj.setMesh(mesh)
            obj.mesh = mesh;
        end
        
        function setSurface(obj, surface)
            % setSurface assigns a surface onto the reflector
            %   obj.setSurface(surface)
            obj.surface = surface;
        end
        
        function [x,y,z] = getMeshPointsI(obj,option)
            % getMeshPointsI returns the coordinates of the integration
            %   points of the mesh, in the reflector CS. The reflector must
            %   have been meshed with meshReflector first.
            %   [x,y,z] = obj.getMeshPointsI()
            %   [x,y,z] = obj.getMeshPointsI('local') returns the
            %   coordinates in the local coordinate system
            if nargin == 2 && strcmp(option, 'local')
                [x,y,z] = obj.mesh.getPointsI(); 
                %disp('getting from local I')
            elseif nargin == 2 && strcmp(option, 'global')
                [x,y,z] = obj.mesh.getPointsI();
                CS = CSlocal(0,0,0,[1,0,0],[0,1,0],[0,0,1]);
                [x,y,z] = CS.pos2CS(x,y,z,obj);
                disp('getting from global I')
            elseif nargin == 1
                %warning('to avoid confusion, please specify global or local')
                [x,y,z] = obj.mesh.getPointsI();
                %disp('getting from local I')
            else
                error('get mesh points in which CS???')
            end
        end
        
        function [x1,x2,J] = getGridPoints(obj)
            % getGridPoints only works with cylindrical and cartesian
            % meshes. It returns the integration grid and the jacobian at
            % those points
            [x1,x2,J] = obj.mesh.getGridPoints();
        end
        
        function [x,y,z] = getMeshPointsE(obj,option)
            % getMeshPointsI returns the coordinates of the vertices
            %   of the mesh, in the reflector CS. The reflector must
            %   have been meshed with meshReflector first.
            %   [x,y,z] = obj.getMeshPointsE()
            %   [x,y,z] = obj.getMeshPointsE('local') returns the
            %   coordinates in the local coordinate system
            if nargin == 2 && strcmp(option, 'local')
                [x,y,z] = obj.mesh.getPointsE();   
                %disp('getting from local E')
            elseif nargin == 2 && strcmp(option, 'global')
                [x,y,z] = obj.mesh.getPointsE();
                CS = CSlocal(0,0,0,[1,0,0],[0,1,0],[0,0,1]);
                [x,y,z] = CS.pos2CS(x,y,z,obj);
                disp('getting from global E')
            elseif nargin == 1
                warning('to avoid confusion, please specify global or local')
                [x,y,z] = obj.mesh.getPointsE();
                %disp('getting from local E')
            else
                error('get mesh points in which CS???')
            end
        end
        
        function [n_x, n_y, n_z, N] = getMeshNormals(obj, option)
            % getMeshNormals returns the normals and Jacobian(?) at the
            %   integration points of the mesh. The reflector must have been
            %   meshed with meshReflector and given 3D shape with
            %   depthReflector first. The normals are expressed in the
            %   reflector CS
            %   [n_x, n_y, n_z, N] = obj.getMeshNormals()
            %   [n_x, n_y, n_z, N] = obj.getMeshNormals('local') returns
            %   the normals in the CS of the reflector
            [n_x, n_y, n_z, N] = obj.mesh.getNormals();
            if nargin == 1
            end
            if nargin == 2 && strcmp(option, 'local') 
                 %disp('getting normals from local')
            elseif nargin == 2 && strcmp(option, 'global') 
                CS = CSlocal(0,0,0,[1,0,0],[0,1,0],[0,0,1]);
                [n_x, n_y, n_z] = CS.vec2CS(n_x, n_y, n_z,obj);
                disp('getting normals from global')
            elseif nargin == 1
                %warning('to avoid confusion, please specify global or local')
                %disp('getting normals from local')
            else
                error('get mesh points in which CS???')
            end
        end
        
        function depthReflector(obj)
            % depthReflector calculates the z position and normals of the
            %   reflector. The reflector must have been meshed with 
            %   meshReflector first, and given a surface with setSurface. 
            if ~obj.mesh.ismeshed || isempty(obj.surface)
                warning('Reflector must be meshed and given a surface with meshReflector before applying profile')
                return
            end

            % Retrieve (x,y) coordinates from mesh
            [xE,yE,~] = obj.mesh.getPointsE();

            % Compute z coordinate
            zE = obj.surface.getSurf(xE, yE);

            % Set z coordinate
            obj.mesh.setZE(zE); 

            % Retrieve (x,y) coordinates from mesh
            [xI,yI,~] = obj.mesh.getPointsI();

            % Compute z coordinate
            zI = obj.surface.getSurf(xI, yI);

            % Set z coordinate
            obj.mesh.setZI(zI);

            % Compute normals and N
            [n_x, n_y, n_z, N] = obj.surface.getNorm(xI, yI);

            % Set normals and N
            obj.mesh.setNormals(n_x, n_y, n_z, N);
        end
        
        function meshReflector(obj)
            % meshReflector meshes the reflector. Must have been given a
            % mesh with setMesh first.
            if isempty(obj.mesh)
                warning('Reflector must given a mesh type with setMesh before meshing')
                return
            end
            obj.mesh.meshNow(obj.rim,obj.surface.getLines());
        end
                
        function taper_angle = getTaperAngle(obj, feed,low_up)

            % Retrieve feed position in global CS
            position_feed=feed.position;

            % Retrieve rim center in global CS
            [x0, y0, z0] = obj.getRimCenter('global');

            % Select upper or lower edge of the reflector
            if strcmp(low_up,'low')
                x1 = obj.rim.x0-obj.rim.a1; 
                y1 = obj.rim.y0;
            elseif strcmp(low_up,'up')
                x1 = obj.rim.x0+obj.rim.a1; 
                y1 = obj.rim.y0;
            end

            % Compute z coordinate
            z1 = obj.surface.getSurf(x1,y1);

            % Define global CS
            CSglobal  = CSlocal(0,0,0,[1 0 0],[0 1 0],[0 0 1]);

            % Convert from reflector CS to gloal CS
            [x1,y1,z1] = CSglobal.pos2CS(x1,y1,z1,obj);

            % Compute angle
            D = [x0,y0,z0]';        % Reflector center (global CS)
            C = [x1,y1,z1]';        % Reflector edge (global CS)
            A = position_feed;      % Feed center (global CS)
            taper_angle = acos(dot(D-A,C-A)/norm(D-A)/norm(C-A));
        end
        
        function disp2fig(obj,ax,CSref)
            % This method can be invoked to plot the reflector
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
                elseif isa(ax,'matlab.ui.control.UIAxes')
                    hold(ax,'on');
                else
                    figure(ax); hold on;
                end
            end
            [xp, yp, zp] = obj.getMeshPointsE('local');
            xp = obj.nanSurfValue(xp);
            yp = obj.nanSurfValue(yp);
            zp = obj.nanSurfValue(zp);
            [xp,yp,zp] = CSref.pos2CS(xp,yp,zp,obj);
            if size(xp,2) >1
                if isa(ax,'matlab.ui.control.UIAxes')
                    surf(ax, xp, yp, zp, 0*zp);
                else
                    surf(xp, yp, zp, 0*zp);
                end
            else
                if isa(ax,'matlab.ui.control.UIAxes')
                    trisurf(ax, obj.mesh.tri, xp,yp,zp,0*zp);
                else
                    trisurf(obj.mesh.tri, xp,yp,zp,0*zp);
                end
            end
            if isa(ax,'matlab.ui.control.UIAxes')
                hold(ax,'on'); axis(ax,'equal');
            else
                hold on; axis equal;
            end
            obj.plotTriade(ax, CSref);
            
            if isfield(obj.plotOptions, 'nplot') && ~isnan(obj.plotOptions.nplot)
                [xi, yi, zi] = obj.getMeshPointsI('local');
                xi = obj.nanSurfValue(xi(:));
                yi = obj.nanSurfValue(yi(:));
                zi = obj.nanSurfValue(zi(:));
                [xi,yi,zi] = CSref.pos2CS(xi,yi,zi,obj);
                toPlotId = randperm(length(xi),obj.plotOptions.nplot);
                [n_x, n_y, n_z, ~] = obj.getMeshNormals();
                [n_x, n_y, n_z] = CSref.vec2CS(n_x(:), n_y(:), n_z(:),obj);
                quiver3(xi(toPlotId),yi(toPlotId),zi(toPlotId),...
                    n_x(toPlotId),n_y(toPlotId),n_z(toPlotId),obj.CSscale/2);
            end
        end
        
    end
    
end
