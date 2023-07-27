classdef Currents2Farfield<Sources2Farfield
    % Currents2Farfield allows to compute the farfield EM field generated
    % by the equivalent surface currents on the reflector on the angular
    % points defined in the coverage object.
    % 
    % Authors : Salvador Mercader Pellicer, Louis Dufour
    % Revisions : v0.1.0 - 25/10/2018

    properties
        
    end

    methods

        function obj = Currents2Farfield()
            % Constructor
            
        end            
        
        function calculateFarField(obj,surfCurrents,coverage,copolar_component)

            % Call function to compute FF
            [obj.D_V_dB, obj.D_H_dB,obj.D_RHCP_dB, obj.D_LHCP_dB,...
                obj.E_Theta, obj.E_Phi, obj.E_V, obj.E_H, obj.E_RHCP, obj.E_LHCP] = ...
                integrateFFutility(surfCurrents, coverage);

            % Assign surfCurrents and coverage to object
            obj.sources=surfCurrents;
            obj.coverage=coverage;
            obj.copolar_component=copolar_component;

        end
    end
end

%% Auxiliary functions

function [D_V_dB, D_H_dB,D_RHCP_dB, D_LHCP_dB,...
    E_Theta, E_Phi, E_V, E_H, E_RHCP, E_LHCP] =  integrateFFutility(currents, coverage)
    % This function solves the FF integral
    
    % Store variables
    reflector = currents.reflector;
    type = reflector.mesh.intType;  % mesh type 'cartesian','cylindrical', etc..
    Prad = currents.Prad;           % Prad of the feed
    
    % For the mesh elements outside the rim of the reflector the surface
    % currents are set to NaN, so not_nan_mesh corresponds to the mask
    % containing the mesh points inside the rim.
    not_nan_mesh=(~isnan(currents.J_x))&(~isnan(currents.J_y))&(~isnan(currents.J_z));
    
    % Remove NaN elements
    J_x = currents.J_x(not_nan_mesh);
    J_y = currents.J_y(not_nan_mesh);
    J_z = currents.J_z(not_nan_mesh);
    
    % Store k_0, eta_0 and phase_yn
    k_0 = currents.k_0;
    eta_0 = currents.eta_0;
    phase_yn = currents.phase_yn;   % 1 if phase added when creating feed; 0 if phase added before performing the integral (0 always when using Ludwig Method)
    
    % if phase_yn==1 the surface currents already include the spherical phase 
    % term; if phase_yn==0 we must include the spherical phase term, so we need
    % the distance between the feed and each mesh points, which is stored in
    % rho_s.
    if phase_yn==1
        rho_s=0;
    else
        rho_s = currents.rho_s(not_nan_mesh);         % Distance of each point in the integration grid with respect to the feed position
    end
    clear currents;
    
    % Extract farfield points
    u = coverage.U;                     % u grid
    v = coverage.V;                     % v grid
    w = coverage.W;                     % w=sqrt(1-(u.^2+v.^2))=cos(theta)
    sinPhi_uv = coverage.sinPhi;        % sin(phi)
    sinTheta_uv = coverage.sinTheta;    % sin(theta)
    cosPhi_uv = coverage.cosPhi;        % cos(phi)
    cosTheta_uv=coverage.cosTheta;      % cos(theta)
    
    % Extract reflector mesh
    [xI, yI, zI] = reflector.getMeshPointsI('local');
    [~,~,~,N] = reflector.getMeshNormals('local');      % N=sqrt(nx^2+ny^2+nz^2)=sqrt(1+d_x(f)^2+d_y(f)^2) is the norm of the normal vector to the surface z=f(x,y)
    
    % Remove points outside the rim
    xI=xI(not_nan_mesh);
    yI=yI(not_nan_mesh);
    zI=zI(not_nan_mesh);
    N=N(not_nan_mesh);
    
    % Extract integration parameters: weights and Jacobian
    switch type
        case {'cartesian','cylindrical'}
            % rxG,phiyG, rx_single and phiy_single are used in the Ludwig
            % integration. weights and Jacobian are used in the trapz
            % integration.
            weights = reflector.mesh.weights(not_nan_mesh);
            [~,~,Jacobian] = reflector.getGridPoints();   % Get the x,y or r,phi points on the integration mesh and Jacobian
            if numel(Jacobian)>1
                Jacobian=Jacobian(not_nan_mesh);
            end
        case 'triGauss'
            weights = reflector.mesh.weights(:);
            Jacobian = 1;
        case 'PolarGRASP'
            weights = reflector.mesh.weights(not_nan_mesh);
            [~,~,Jacobian] = reflector.getGridPoints();
            Jacobian=Jacobian(not_nan_mesh);
    end
    
    % Extract FF points grid
    [N_v,N_u]=size(u);
    
    % Initialize some auxiliary variables 
    FTx=nan(N_v,N_u);
    FTy=nan(N_v,N_u);
    FTz=nan(N_v,N_u);
    
    % Avoid computation of NaN values in the u-v plane
    not_NaN_mask=~isnan(u);
    u=u(not_NaN_mask);
    v=v(not_NaN_mask);
    w=w(not_NaN_mask);
    N_not_NaN=length(u);
    FTx_aux=nan(N_not_NaN,1);
    FTy_aux=nan(N_not_NaN,1);
    FTz_aux=nan(N_not_NaN,1);
    
    % N is the norm of the n vector obtained as n=cross(v1,v2) where v1
    % and v2 are the derivatives of a surface along the two coordinates
    % defining it (see surface integral chapter in a textbook)
    J_x_N=J_x.*N.*exp(-1i*k_0*(rho_s)*(1-phase_yn));
    J_y_N=J_y.*N.*exp(-1i*k_0*(rho_s)*(1-phase_yn));
    J_z_N=J_z.*N.*exp(-1i*k_0*(rho_s)*(1-phase_yn));
    switch type
        case {'cartesian','cylindrical','triGauss','PolarGRASP'}
            parfor n=1:N_not_NaN

                Q=exp(1i*k_0*(u(n)*xI+v(n)*yI+w(n)*zI));

                funx=J_x_N.*Q;
                funy=J_y_N.*Q;
                funz=J_z_N.*Q;

                FTx_aux(n) = weights.'*(funx(:).*Jacobian);
                FTy_aux(n) = weights.'*(funy(:).*Jacobian);
                FTz_aux(n) = weights.'*(funz(:).*Jacobian);
            end
    end
    
    FTx(not_NaN_mask)=FTx_aux;
    FTy(not_NaN_mask)=FTy_aux;
    FTz(not_NaN_mask)=FTz_aux;
    clear FTx_aux FTy_aux FTz_aux;
    
    %% convert to coverage CS
    [x_reflector,y_reflector,z_reflector]=coverage.pos2CS(0,0,0,reflector);
    reflector_to_coverage_CS_converter=exp(1i*k_0*(x_reflector*coverage.U+y_reflector*coverage.V+z_reflector*coverage.W));
    FTx=FTx.*reflector_to_coverage_CS_converter;
    FTy=FTy.*reflector_to_coverage_CS_converter;
    FTz=FTz.*reflector_to_coverage_CS_converter;
    
    %% Field computation
    % Compute electric field
    E_Theta = -1i*k_0*eta_0/(4*pi).*(cosTheta_uv.*cosPhi_uv.*FTx+cosTheta_uv.*sinPhi_uv.*FTy-sinTheta_uv.*FTz);
    E_Phi   = -1i*k_0*eta_0/(4*pi).*(-sinPhi_uv.*FTx+cosPhi_uv.*FTy);
    E_V     = E_Theta.*cosPhi_uv-E_Phi.*sinPhi_uv;
    E_H     = E_Theta.*sinPhi_uv+E_Phi.*cosPhi_uv;
    E_RHCP  = 1/sqrt(2)*(E_V+1i*(E_H));
    E_LHCP  = 1/sqrt(2)*(E_V-1i*(E_H));
    
    % Compute gain
    D_V_dB    = 10*log10(abs(E_V)   .^2/(2*eta_0)*4*pi/Prad);
    D_H_dB    = 10*log10(abs(E_H)   .^2/(2*eta_0)*4*pi/Prad);
    D_RHCP_dB = 10*log10(abs(E_RHCP).^2/(2*eta_0)*4*pi/Prad);
    D_LHCP_dB = 10*log10(abs(E_LHCP).^2/(2*eta_0)*4*pi/Prad);
    
    % Apply distance factor
    if ~isempty(coverage.r)
        distance_factor  = exp(-1i*k_0*coverage.r)./coverage.r;
        
        E_Theta = E_Theta .*distance_factor;
        E_Phi   = E_Phi   .*distance_factor;
        E_V     = E_V     .*distance_factor;
        E_H     = E_H     .*distance_factor;
        E_RHCP  = E_RHCP  .*distance_factor;
        E_LHCP  = E_LHCP  .*distance_factor;
    end

end