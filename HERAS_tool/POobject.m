classdef (Abstract) POobject < matlab.mixin.Copyable
    % POobject is the abstract class for any geometrical object
    %   It is the mist basic element class, where each object has a
    %   position and a direction, as well as scale for CS plotting
    %   It provides basic functions, such as coordinate system
    %   transformations
    %   Author : Louis Dufour                                             
	%   Revisions : 0.1.0 - 25/10/2018  
    %
    % POobject properties:
    %   position    - position of the object
    %   direction   - direction of the object (orhtonormal coordinate
    %                 system)
    %   CSscale     - scale used to plot the coordinate system of the
    %                 object  
    %
    % POobject methods:
    %   pos2CS      - convert points from another object's coordinate
    %                 system into the current object cordinate system
    %   vec2CS      - converts vectors from another CS to its own CS
    %   getCenter   - returns the center x,y,z of the object
    %   getCStriade - return the local CS, expressed in the global CS
    properties
        position = nan(3,1);  % 3 by 1 vector, which represents the position of the object in 3D
        direction = nan(3,3); % 3 by 3 matrix made of the three orhtonormal coordinate system axis 
                              % [ex ey ez], with ei a 3 by 1 vector
        CSscale = 0.1;        % scale used for plotting the coordinate system of the object
        CSind;
        plotOptions;
    end
    
    methods (Abstract)
        disp2fig(obj,ax,CSref);
    end
    
    methods
        function [coord_1x, coord_1y, coord_1z] = pos2CS(obj,coord_2x, coord_2y, coord_2z, obj2)
            % pos2CS converts the coordinates from the coordinates system of
            %   another object obj2 into the coordinate system of the current
            %   object.
            %   [x1, y1, z1] = obj1.pos2CS(x2, y2, z2, obj2) returns the
            %   coordinates, in the coordinate system of obj1, of the
            %   points given by [x2,y2,z2] in the CS of obj2
            %   Ex : we have a feed and a reflector, both instances of the
            %   POobject class. The reflector has, in its own coodinate
            %   system, the x2,y2,z2 integration points. To calculate the
            %   surface current at those points, we need to express them in
            %   the feed coordinate system:
            %   [x1, y1, z1] = feed.pos2CS(x2, y2, z2, reflector)
            %
            % See POReflector, POFeed 

            %% Francesco Lisi update - 01/03/23
            
%             % Compute total number of points
%             N = numel(coord_2x);

            % definition of a structure with obj CS position and direction
            CS1.position=obj.position;
            CS1.direction=obj.direction;

            % definition of a structure with obj2 CS position and direction
            CS2.position=obj2.position;
            CS2.direction=obj2.direction;

            % Transformation
            [coord_1x, coord_1y, coord_1z] = pos2CSfcn(CS2,CS1,coord_2x, coord_2y, coord_2z);


%             % Extract CS origin
%             pos1 = obj.position;
%             pos2 = obj2.position;
% 
%             % Extract CS direction matrix
%             dir1 = obj.direction;
%             dir2 = obj2.direction;
%             
%             % Extract input coordinate size
%             coord_size = size(coord_2x);
% 
%             % creation of a matrix whose column correspond to the three
%             % coordinates of each point in the obj2 CS
%             points_2=[coord_2x(:).';coord_2y(:).';coord_2z(:).'];
% %             points_2=zeros(3,N);
% %             points_2(1,:) = reshape(coord_2x,[1,N]);
% %             points_2(2,:) = reshape(coord_2y,[1,N]);
% %             points_2(3,:) = reshape(coord_2z,[1,N]);
% 
%             % Conversion of the points from obj2 CS to obj CS
%             points_1 = dir1'*(dir2*points_2+pos2-pos1);
% 
%             % Output extraction
%             coord_1x = reshape(points_1(1,:),coord_size);
%             coord_1y = reshape(points_1(2,:),coord_size);
%             coord_1z = reshape(points_1(3,:),coord_size);
% 

%          	pos1 = obj.position;
%             pos2 = obj2.position;
%             dir1 = obj.direction;
%             dir2 = obj2.direction;
%             %R = mat1(:,1) *  mat2(:,1)' +  mat1(:,2) *  mat2(:,2)' +  mat1(:,3) *  mat2(:,3)';
%             R = inv(dir1)*dir2;
%             coord_1x = nan(size(coord_2x));
%             coord_1y = nan(size(coord_2x));
%             coord_1z = nan(size(coord_2x));
%             dpos = dir2\(pos2-pos1);
%             
%             aux1=[coord_2x(:), coord_2y(:), coord_2z(:)].'+dpos*ones(size(dpos,2),length(coord_2x(:)));
%             aux2=R*aux1;
%             aux3=reshape(aux2,3,size(coord_2x,1),size(coord_2x,2));
%             coord_1x=reshape(squeeze(aux3(1,:,:)),size(coord_2x,1),size(coord_2x,2));
%             coord_1y=reshape(squeeze(aux3(2,:,:)),size(coord_2x,1),size(coord_2x,2));
%             coord_1z=reshape(squeeze(aux3(3,:,:)),size(coord_2x,1),size(coord_2x,2));
            
%             tic;
%             for ii = 1:size(coord_2x,1)
%                 for jj = 1:size(coord_2x,2)
%                     coord_1t = R*([coord_2x(ii,jj); coord_2y(ii,jj); coord_2z(ii,jj)] + dpos);
%                     coord_1x(ii,jj) = coord_1t(1);
%                     coord_1y(ii,jj) = coord_1t(2);
%                     coord_1z(ii,jj) = coord_1t(3);
%                 end
%             end
%             toc;
            
        end

        function [coord_1x, coord_1y, coord_1z] = vec2CS(obj,coord_2x, coord_2y, coord_2z, obj2)
            % vec2CS converts the vectors from the coordinates system of
            %   another object obj2 into the coordinate system of the current
            %   object.
            %   Used for instance to calculate the H fields from the feed
            %   illumination of the reflector
            %{
            %   Can also accept two parameters (in which case the first one
            %   is a N x 3 matrix and the second obj2), and output only one
            %   argument (a Nx3 matrix)
            if nargin == 3
                obj2 = coord_2y;
                coord_2y = coord_2x(:,2);
                coord_2z = coord_2x(:,3);
                coord_2x = coord_2x(:,1);
            end
            %}

            %% Francesco Lisi update - 01/03/23
            
            % definition of a structure with obj CS position and direction
            CS1.position=zeros(3,1);
            CS1.direction=obj.direction;

            % definition of a structure with obj2 CS position and direction
            CS2.position=zeros(3,1);
            CS2.direction=obj2.direction;

            % Transformation
            [coord_1x, coord_1y, coord_1z] = pos2CSfcn(CS2,CS1,coord_2x, coord_2y, coord_2z);
            
%             % Extract CS direction matrix
%             dir1 = obj.direction;
%             dir2 = obj2.direction;
%             
%             % Extract input vector size
%             coord_size = size(coord_2x);
% 
%             % creation of a matrix whose column correspond to the three
%             % components of each point in the obj2 CS
%             points_2=[coord_2x(:).';coord_2y(:).';coord_2z(:).'];
% %             points_2=zeros(3,N);
% %             points_2(1,:) = reshape(coord_2x,[1,N]);
% %             points_2(2,:) = reshape(coord_2y,[1,N]);
% %             points_2(3,:) = reshape(coord_2z,[1,N]);
% 
%             % Conversion of the points from obj2 CS to obj CS
%             points_1 = dir1'*(dir2*points_2);
% 
%             % Output extraction
%             coord_1x = reshape(points_1(1,:),coord_size);
%             coord_1y = reshape(points_1(2,:),coord_size);
%             coord_1z = reshape(points_1(3,:),coord_size);
            
%             mat2 = obj2.direction;
%             mat1 = obj.direction;
%             %R = mat1(:,1) *  mat2(:,1)' +  mat1(:,2) *  mat2(:,2)' +  mat1(:,3) *  mat2(:,3)'
%             R = inv(mat1)*mat2;
%             coord_1x = nan(size(coord_2x));
%             coord_1y = nan(size(coord_2x));
%             coord_1z = nan(size(coord_2x));
%             dpos = 0;
%             
%             aux1=[coord_2x(:), coord_2y(:), coord_2z(:)].'+dpos*ones(size(dpos,2),length(coord_2x(:)));
%             aux2=R*aux1;
%             aux3=reshape(aux2,3,size(coord_2x,1),size(coord_2x,2));
%             coord_1x=reshape(squeeze(aux3(1,:,:)),size(coord_2x,1),size(coord_2x,2));
%             coord_1y=reshape(squeeze(aux3(2,:,:)),size(coord_2x,1),size(coord_2x,2));
%             coord_1z=reshape(squeeze(aux3(3,:,:)),size(coord_2x,1),size(coord_2x,2));

            
%             for ii = 1:size(coord_2x,1)
%                 for jj = 1:size(coord_2x,2)
%                     coord_1t = R*([coord_2x(ii,jj); coord_2y(ii,jj); coord_2z(ii,jj)] +dpos);
%                     coord_1x(ii,jj) = coord_1t(1);
%                     coord_1y(ii,jj) = coord_1t(2);
%                     coord_1z(ii,jj) = coord_1t(3);
%                 end
%             end
            %{
            if nargout == 1
                coord_1x = [coord_1x coord_1y coord_1z];
            end
            %}
        end
        
        function [x,y,z] = getCenter(obj)
            % getCenter returns the object's center coordinates (in the
            % global CS)
            %   [x,y,z] = obj.getCenter()
            x = obj.position(1);
            y = obj.position(2);
            z = obj.position(3);
        end
        
        function [ex, ey, ez] = getCStriade(obj)
            % getCStriade returns the object's orthonormal coordinates system
            % (in the global CS)
            %   [ex, ey, ez] = obj.getCStriade()
            ex = obj.direction(:,1);
            ey = obj.direction(:,2);
            ez = obj.direction(:,3);
        end
        
        function rotateLocal(obj,dir,angle)
            % ROTATELOCAL rotates the object in its local CS, along on of
            % its axis (DIR = 'x', or 'y', or 'z') of an angle ANGLE (in
            % radians)
            switch dir
                case 'x'
                    u = obj.direction(:,1);
                case 'y'
                    u = obj.direction(:,2);
                case 'z'
                    u = obj.direction(:,3);
            end
            x = u(1); y  = u(2); z = u(3);
            M = [0 -z y;z 0 -x;-y x 0]*sin(angle)+(eye(3)-u*u')*cos(angle)+u*u';
            obj.direction = M*obj.direction;
        end
        
        function displaceLocal(obj,displ)
            % DISPLACELOCAL displace the object in its local CS along on of
            % its axis of a distance DISPL (DISPL=[x_displ;y_displ;z_displ]
            obj.position=obj.position+obj.direction*displ;
        end
        
        function plotTriade(obj,ax, CSref)
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
            
            %[x,y,z] = obj.getCenter();
            %[ex,ey,ez] = obj.getCStriade();
            [x,y,z] = CSref.pos2CS(0,0,0,obj);
            [exx, exy, exz] = CSref.vec2CS(1,0,0, obj);
            [eyx, eyy, eyz] = CSref.vec2CS(0,1,0, obj);
            [ezx, ezy, ezz] = CSref.vec2CS(0,0,1, obj);
            ex = [exx; exy; exz]; ey = [eyx; eyy; eyz]; ez = [ezx; ezy; ezz];
            ex = ex*obj.CSscale; ey = ey*obj.CSscale; ez = ez*obj.CSscale;
            
            if isa(ax,'matlab.ui.control.UIAxes')
                if isempty(obj.CSind)
                    text(ax,x+ex(1),y+ex(2),z+ex(3),'x');
                    text(ax,x+ey(1),y+ey(2),z+ey(3),'y');
                    text(ax,x+ez(1),y+ez(2),z+ez(3),'z');
                else
                    text(ax,x+ex(1),y+ex(2),z+ex(3),['x_{' obj.CSind '}']);
                    text(ax,x+ey(1),y+ey(2),z+ey(3),['y_{' obj.CSind '}']);
                    text(ax,x+ez(1),y+ez(2),z+ez(3),['z_{' obj.CSind '}']);
                end
                quiver3(ax,x,y,z,ex(1),ex(2),ex(3),'r');
                quiver3(ax,x,y,z,ey(1),ey(2),ey(3),'g');
                quiver3(ax,x,y,z,ez(1),ez(2),ez(3),'b');
            else
                if isempty(obj.CSind)
                    text(x+ex(1),y+ex(2),z+ex(3),'x');
                    text(x+ey(1),y+ey(2),z+ey(3),'y');
                    text(x+ez(1),y+ez(2),z+ez(3),'z');
                else
                    text(x+ex(1),y+ex(2),z+ex(3),['x_{' obj.CSind '}']);
                    text(x+ey(1),y+ey(2),z+ey(3),['y_{' obj.CSind '}']);
                    text(x+ez(1),y+ez(2),z+ez(3),['z_{' obj.CSind '}']);
                end
                quiver3(x,y,z,ex(1),ex(2),ex(3),'r');
                quiver3(x,y,z,ey(1),ey(2),ey(3),'g');
                quiver3(x,y,z,ez(1),ez(2),ez(3),'b');
            end
            
            
        end
    end
end