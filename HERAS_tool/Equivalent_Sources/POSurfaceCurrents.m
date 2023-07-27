classdef POSurfaceCurrents < FarfieldEquivalentSources
    % POSurfaceCurrents is used to associate a reflector with a feed 
    % calculate the surface currents create by the feed or
    % the subreflector on the reflector. The first parameter is a
    % reflector, and the second parameter is a feed
    % To use it, first create the object with POSurfaceCurrents(reflector,
    % feed), and then calculate the currents by running
    % calculateCurrents().
    % Authors : Salvador Mercader Pellicer, Louis Dufour
    % Revisions : v0.2.0 - 25/10/2018
    % Revisions : v0.3.0 - 13/02/2023 - Author: Francesco Lisi
    %             - POSurfaceCurrents is the child class of the class
    %               FarfieldEquivalentSources;
    %             - Fixed a bug in the computation of the reflected fields.


    properties
        J_x, J_y, J_z;
    end

    methods

        function obj = POSurfaceCurrents(reflector, feed, options)
            % Constructor
            obj.reflector = reflector;
            obj.feed = feed;
            obj.k_0 = feed.k_0;
            obj.phase_yn = feed.phase_yn;
            obj.eta_0 = feed.eta_0;
            obj.Prad = feed.Prad;
            if nargin == 3
                obj.options  = options;
            else
                obj.options = [];
            end

        end
                
        
        function calculateCurrents(obj)
            % This method computes the surface currents on the reflector
            % mesh

            % Get the integration points
            [x,y,z] = obj.reflector.getMeshPointsI('local');

            % Compute mask with points on the mesh that are outside the
            % reflector rim.
            outside_mask=~obj.reflector.InsideRimMask();

            % Set values outside the rim to NaN
            x(outside_mask)=NaN;
            y(outside_mask)=NaN;
            z(outside_mask)=NaN;

            % store mesh size
            mesh_size=size(x);

            % Express mesh points in the feed coordinate system
            [x_s,y_s,z_s] = obj.feed.pos2CS(x,y,z,obj.reflector);
            clear x y z;
            
            % Compute mesh point distance from feed
            [~,~,obj.rho_s]=cart2sph(x_s,y_s,z_s);

            % STEP ONE: GET THE INCIDENT FIELDS
            % Compute incident fields in the feed CS
            [E_x_s,E_y_s,E_z_s,H_x_s,H_y_s,H_z_s] = obj.feed.illuminate(x_s,y_s,z_s);
            clear x_s y_s z_s;
            
            N=1;

            % Convert the vectors to the reflector CS
            [obj.H_x_i,obj.H_y_i,obj.H_z_i] = obj.reflector.vec2CS(H_x_s,H_y_s,H_z_s, obj.feed);
            clear H_x_s H_y_s H_z_s;
            [obj.E_x_i,obj.E_y_i,obj.E_z_i] = obj.reflector.vec2CS(E_x_s,E_y_s,E_z_s, obj.feed);
            clear E_x_s E_y_s E_z_s;
            
            % STEP TWO: GET THE REFLECTED FIELDS
            % Create incident field arrays
            E_i=permute(reshape([obj.E_x_i(:).';obj.E_y_i(:).';obj.E_z_i(:).'],[3,N,mesh_size]),[1,3,4,2]);
            H_i=permute(reshape([obj.H_x_i(:).';obj.H_y_i(:).';obj.H_z_i(:).'],[3,N,mesh_size]),[1,3,4,2]);

            % Get the normals at the integration points (in the reflector CS)
            [n_x, n_y, n_z, ~] = obj.reflector.getMeshNormals('local');

            n=reshape([n_x(:).';n_y(:).';n_z(:).'],[3,mesh_size]);
            clear n_x n_y n_z;

            % Computation of the reflected fields
            E_r=-E_i+2*(sum(E_i.*n,1)).*n;
            H_r=+H_i-2*(sum(H_i.*n,1)).*n;

            % Store reflected fields
            obj.E_x_r = permute(reshape(E_r(1,:),[mesh_size,N]),[3,1,2]);
            obj.E_y_r = permute(reshape(E_r(2,:),[mesh_size,N]),[3,1,2]);
            obj.E_z_r = permute(reshape(E_r(3,:),[mesh_size,N]),[3,1,2]);
            obj.H_x_r = permute(reshape(H_r(1,:),[mesh_size,N]),[3,1,2]);
            obj.H_y_r = permute(reshape(H_r(2,:),[mesh_size,N]),[3,1,2]);
            obj.H_z_r = permute(reshape(H_r(3,:),[mesh_size,N]),[3,1,2]);
            
            % STEP3 : GET THE SURFACE CURRENTS
            % Calculate the suface currents at the integration points
            J=2*cross(repmat(n,[1,1,1,N]),H_i,1);

            % Store currents
            obj.J_x=permute(reshape(J(1,:),[mesh_size,N]),[3,1,2]);
            obj.J_y=permute(reshape(J(2,:),[mesh_size,N]),[3,1,2]);
            obj.J_z=permute(reshape(J(3,:),[mesh_size,N]),[3,1,2]);
            clear J E_r H_r E_i H_i;

            % Remove singleton dimension
            obj.E_x_r=squeeze(permute(obj.E_x_r,[2,3,1]));
            obj.E_y_r=squeeze(permute(obj.E_y_r,[2,3,1]));
            obj.E_z_r=squeeze(permute(obj.E_z_r,[2,3,1]));
            obj.H_x_r=squeeze(permute(obj.H_x_r,[2,3,1]));
            obj.H_y_r=squeeze(permute(obj.H_y_r,[2,3,1]));
            obj.H_z_r=squeeze(permute(obj.H_z_r,[2,3,1]));
            obj.J_x=squeeze(permute(obj.J_x,[2,3,1]));
            obj.J_y=squeeze(permute(obj.J_y,[2,3,1]));
            obj.J_z=squeeze(permute(obj.J_z,[2,3,1]));

        end

    end

end