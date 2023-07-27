classdef IdealGaussianSourceFeed<POFeedSingle
% Author : Francesco Lisi
% Revisions : v0.0.1 - 30/06/2023
% -------------------------------------------------------------------------
% This function is a completely redesigned update of the class 
% HuygensSourceFeed (by Salvador Mercader Pellicer, Louis Dufur) to 
% generate a gaussian source in th NF. The radiated power of a gaussian 
% source can be computed in closed form by solving the radiation integral. 
% As a consequence, the numerical integration to compute the power has been 
% removed.
% -------------------------------------------------------------------------
% PROPERTIES:
%   tapering_dB     : tapering value exressed in dB.
%   tapering_angle  : tapering angle expressed in rad. This corresponds to
%                     the theta angle where the FF amplitude drops of 
%                     tapering_dB with respect to the broadside direction.
%   LPCP            : 'LP' for linear and 'CP' for circular polarization.
%   LHRH            : 'LH' for left-hand and 'RH' for right-hand.
%   orientation     : anticlockwise rotation with respect to x axis in rad.
%   r,theta,phi     : polar coordinates where the field is evaluated
%   phase_yn        : 'True'  -> spherical phase term included;
%                     'False' -> spherical phase term not included;
%   k_0             : free space wavenumber
%   eta_0           : free space impedance
%   Prad            : radiated power
%   forceFarField   : if 'True' the FF expression of a gaussian source is
%                     used to compute the EM field


    properties
        tapering_dB,tapering_angle,LPCP,LHRH,orientation;
        forceFarField = true;
    end

    methods

        function obj = IdealGaussianSourceFeed(x,y,z,ex,ey,ez)
            % Constructor
            % (x,y,z) is the position of the feed CS origin expressed in
            % the global CS.
            % (ex,ey,ez) are the (x,y,z) versor coordinates expressed in
            % the global CS.

            % Normalize the versors
            ex = ex/norm(ex); 
            ey = ey/norm(ey); 
            ez = ez/norm(ez);

            % Check that ex, ey and ez are an orthonormla triade
            if dot(ex,ey)>1e-10 ||dot(ex,ez)>1e-10 ||dot(ey,ez)>1e-10
                warning('The coordinate system defined for the feed is not orthogonal');
                return
            end

            % Set position, direction and default power
            obj.position = [x;y;z];
            obj.direction = [ex(:) ey(:) ez(:)];
            obj.Prad = 4*pi;

        end
        
        function setParameters(obj,tapering_dB,tapering_angle,LPCP,LHRH,orientation,phase_yn,k_0,eta_0)
            % This method sets the fundamental parameter of a Huygens
            % source
            obj.LPCP = LPCP;                    
            obj.orientation = orientation;      
            obj.tapering_dB = tapering_dB;      
            obj.k_0 = k_0;                      
            obj.eta_0 = eta_0;                 
            obj.phase_yn = phase_yn;            
            obj.LHRH = LHRH;          
            obj.tapering_angle = tapering_angle;        
            lambda = 2*pi/obj.k_0;             
            c=physconst('LightSpeed');          
            obj.freq = c/lambda;               
        end
        
        function [E_x_s,E_y_s,E_z_s,H_x_s,H_y_s,H_z_s] = illuminate(obj,x_s,y_s,z_s)
            % This method computes the fields at (x_s,y_s,z_s) expressed in
            % feed CS

            % Variable class
            variable_class = class(x_s);

            % Store mesh size
            mesh_size=size(x_s);
            
            % Compute nan mask
            not_nan_mask=~(isnan(x_s)|isnan(y_s)|isnan(z_s));
            not_nan_mask=find(not_nan_mask==1);

            % Remove NaN elements
            x_s=x_s(not_nan_mask);
            y_s=y_s(not_nan_mask);
            z_s=z_s(not_nan_mask);

            % Initialize fields
            E_x_s=complex(nan(mesh_size,variable_class));
            E_y_s=E_x_s;
            E_z_s=E_x_s;
            H_x_s=E_x_s;
            H_y_s=E_x_s;
            H_z_s=E_x_s;

            % Convert from cartesian to polar coordinates
            [r_s,theta_s,phi_s]=cart2sphEM(x_s,y_s,z_s);
    
            % Compute farfield
            if obj.forceFarField
                
                % Compute fields (FF expression) 
                [E_theta_s,E_phi_s,H_theta_s,H_phi_s]=...
                    GaussianSourceFF(obj.tapering_dB,obj.tapering_angle,obj.LPCP,obj.LHRH,obj.orientation,r_s,theta_s,phi_s,obj.phase_yn,obj.k_0,obj.eta_0,obj.Prad);

                % Convert from cartesian to polar coordinates
                [E_x_s(not_nan_mask),E_y_s(not_nan_mask),E_z_s(not_nan_mask)]=sph2cartEMvec(0,E_theta_s,E_phi_s,theta_s,phi_s);
                clear E_theta_s E_phi_s;
                [H_x_s(not_nan_mask),H_y_s(not_nan_mask),H_z_s(not_nan_mask)]=sph2cartEMvec(0,H_theta_s,H_phi_s,theta_s,phi_s);
                clear H_theta_s H_phi_s;

            else

                % Compute fields (NF expression)
                [E_r_s,E_theta_s,E_phi_s,H_r_s,H_theta_s,H_phi_s]=...
                    GaussianSourceNF(obj.tapering_dB,obj.tapering_angle,obj.LPCP,obj.LHRH,obj.orientation,r_s,theta_s,phi_s,obj.phase_yn,obj.k_0,obj.eta_0,obj.Prad);

                % Convert from cartesian to polar coordinates
                [E_x_s(not_nan_mask),E_y_s(not_nan_mask),E_z_s(not_nan_mask)]=sph2cartEMvec(E_r_s,E_theta_s,E_phi_s,theta_s,phi_s);
                clear E_theta_s E_phi_s E_r_s;
                [H_x_s(not_nan_mask),H_y_s(not_nan_mask),H_z_s(not_nan_mask)]=sph2cartEMvec(H_r_s,H_theta_s,H_phi_s,theta_s,phi_s);
                clear H_theta_s H_phi_s H_r_s;                
            end
            
        end
        
        function setForceFarField(obj, bool)
            obj.forceFarField = bool;
        end
    end
end

%% Auxiliary functions

function [E_r,E_theta,E_phi,H_r,H_theta,H_phi]=GaussianSourceNF(tapering_dB,tapering_angle,LPCP,LHRH,orientation,r,theta,phi,phase_yn,k_0,eta_0,Pt)
% Author : Francesco Lisi
% Revisions : v0.0.1 - 30/06/2023
% -------------------------------------------------------------------------
% This function is a completely redesigned update of the function 
% HuygensSourceImagDispNF (by Salvador Mercader Pellicer) to generate
% a gaussian source in th NF. The radiated power of a gaussian source can 
% be computed in closed form by solving the radiation integral. As a 
% consequence, the numerical integration to compute the power has been 
% removed.
% -------------------------------------------------------------------------
% INPUTS:
%   tapering_dB     : tapering value exressed in dB.
%   tapering_angle  : tapering angle expressed in rad. This corresponds to
%                     the theta angle where the FF amplitude drops of 
%                     tapering_dB with respect to the broadside direction.
%   LPCP            : 'LP' for linear and 'CP' for circular polarization.
%   LHRH            : 'LH' for left-hand and 'RH' for right-hand.
%   orientation     : anticlockwise rotation with respect to x axis in rad.
%   r,theta,phi     : polar coordinates where the field is evaluated
%   phase_yn        : True  -> spherical phase term included;
%                     False -> spherical phase term not included;
%   k_0             : free space wavenumber
%   eta_0           : free space impedance
%   Pt              : radiated power
% -------------------------------------------------------------------------
% OUTPUTS:
%   E_r,E_theta,E_phi : r, theta, and phi components electric field.
%   H_r,H_theta,H_phi : r, theta, and phi components magnetic field.

% FF size
FF_size=size(r);
N_FF=prod(FF_size);

% Set b parameter
b=(20*log10((1+cos(tapering_angle))/2)-tapering_dB)/(20*k_0*(1-cos(tapering_angle))*log10(exp(1)));    

% Create FF point vector
[x,y,z]=sph2cartEM(r(:),theta(:),phi(:));
r_vec=[x,y,z].';
clear x y z;

% The gaussian feed is obtained by displacing a Huygens source of 
% [0;0;-1i*b]
R_vec=r_vec-[0;0;-1i*b];
R=sqrt(sum(R_vec.^2,1));
clear r_vec;

% 90 deg rotation matrix
M_90=[0,-1,0;1,0,0;0,0,1];

% Electric and magnetic dipole orientation
ue=[+cos(orientation);+sin(orientation);0];
um=M_90*ue;

% Compute auxiliary variables
E0  = -1i*sqrt(2*eta_0*Pt)/ComputeI(k_0*b).*exp(-1i*phase_yn*k_0*real(R)).*exp(phase_yn*k_0*(imag(R)-b))./R;
aux1=1i./(k_0*R)+1./(k_0*R).^2;
aux2=1-1i./(k_0*R);

% Compute electric and magnetic field linear gaussian source
E=E0.*( (2*aux1.*sum(R_vec.*ue,1).*R_vec+(1-aux1).*cross(cross(R_vec,repmat(ue,[1,N_FF]),1),R_vec,1))./(R.^2) - aux2.*cross(R_vec,repmat(um,[1,N_FF]),1)./R );
H=E0/eta_0.*( aux2.*cross(R_vec,repmat(ue,[1,N_FF]),1)./R + (2*aux1.*sum(R_vec.*um,1).*R_vec + (1-aux1).*cross(cross(R_vec,repmat(um,[1,N_FF]),1),R_vec,1) )./(R.^2) );

if isequal(LPCP,'CP') 
    
    switch LHRH
        case 'LH'
            m=+1;
        case 'RH'
            m=-1;
    end
        
    % Rotate ue and um of 90 deg
    ue=M_90*ue;
    um=M_90*um;

    % Add +/-i field of the 90 deg rotated gaussian feed to obtain CP
    E=E+1i*m*E0.*( (2*aux1.*sum(R_vec.*ue,1).*R_vec+(1-aux1).*cross(cross(R_vec,repmat(ue,[1,N_FF]),1),R_vec,1))./(R.^2) - aux2.*cross(R_vec,repmat(um,[1,N_FF]),1)./R );
    H=H+1i*m*E0/eta_0.*( aux2.*cross(R_vec,repmat(ue,[1,N_FF]),1)./R + (2*aux1.*sum(R_vec.*um,1).*R_vec + (1-aux1).*cross(cross(R_vec,repmat(um,[1,N_FF]),1),R_vec,1) )./(R.^2) );

    % Divide by sqrt(2) the fields
    E=E/sqrt(2);
    H=H/sqrt(2);
    
end
clear E0 aux 1 aux 2 ue um R R_vec;

% Extract polar component of the fields
[E_r,E_theta,E_phi]=cart2sphEMvec(E(1,:),E(2,:),E(3,:),theta(:).',phi(:).');
[H_r,H_theta,H_phi]=cart2sphEMvec(H(1,:),H(2,:),H(3,:),theta(:).',phi(:).');

% Reshape output
E_r     = reshape(E_r     ,FF_size);
E_theta = reshape(E_theta ,FF_size);
E_phi   = reshape(E_phi   ,FF_size);
H_r     = reshape(H_r     ,FF_size);
H_theta = reshape(H_theta ,FF_size);
H_phi   = reshape(H_phi   ,FF_size);

end

function [E_theta,E_phi,H_theta,H_phi]=GaussianSourceFF(tapering_dB,tapering_angle,LPCP,LHRH,orientation,r,theta,phi,phase_yn,k_0,eta_0,Pt)
% Author : Francesco Lisi
% Revisions : v0.0.1 - 29/06/2023
% -------------------------------------------------------------------------
% This function is a completely redesigned update of the function 
% HuygensSourceImagDisplacement (by Salvador Mercader Pellicer) to generate
% a gaussian source. The radiated power of a gaussian source can e computed
% in closed form by solving the radiation integral. As a consequence, the
% numerical integration to compute the power has been removed.
% -------------------------------------------------------------------------
% INPUTS:
%   tapering_dB     : tapering value exressed in dB.
%   tapering_angle  : tapering angle expressed in rad. This corresponds to
%                     the theta angle where the FF amplitude drops of 
%                     tapering_dB with respect to the broadside direction.
%   LPCP            : 'LP' for linear and 'CP' for circular polarization.
%   LHRH            : 'LH' for left-hand and 'RH' for right-hand.
%   orientation     : anticlockwise rotation with respect to x axis in rad.
%   r,theta,phi     : polar coordinates where the field is evaluated
%   phase_yn        : 'True'  -> spherical phase term included;
%                     'False' -> spherical phase term not included;
%   k_0             : free space wavenumber
%   eta_0           : free space impedance
%   Pt              : radiated power
% -------------------------------------------------------------------------
% OUTPUTS:
%   E_theta,E_phi   : theta and phi components electric field.
%   H_theta,H_phi   : theta and phi components magnetic field.

% Set b parameter
b=(20*log10((1+cos(tapering_angle))/2)-tapering_dB)/(20*k_0*(1-cos(tapering_angle))*log10(exp(1)));    

E0  = -1i*sqrt(2*eta_0*Pt)/ComputeI(k_0*b)*(1+cos(theta)).*exp(k_0*b*(cos(theta)-1)).*exp(-1i*k_0*(r)*phase_yn)./r;

switch LPCP 

    % Linear Polarization
    case 'LP'

        % Compute theta and phi components
        E_theta   = E0.*cos(orientation-phi);
        E_phi     = E0.*sin(orientation-phi);

    % Circular Polarization
    case 'CP' 
        
        % Set m parameter depending on rotation direction (LH or RH)
        switch LHRH
            case 'LH'
                m=+1;
            case 'RH'
                m=-1;
        end
    
        % Update E0
        E0 = E0.*exp(-1i*m*(orientation-phi));
    
        % Compute theta and phi components of the electric field
        E_theta=E0/sqrt(2);
        E_phi=1i*m*E_theta;
    
end
clear E0;

% Spherical components of magnetic field
H_theta   = -E_phi   /eta_0;
H_phi     = +E_theta /eta_0;

end

function I = ComputeI(a)
% This function computes the integral of the function
% (1+cos(theta)).^2.*exp(2*a*cos(theta)).*sin(theta) for theta in [0,pi)
% and phi in [0,2*pi). The exp(2*a) factor is removed and included in the
% E0 expression to avoid inf results
I=sqrt(((2/a)-1/(a^2)+1/(4*a^3))-(1/(4*a^3))*exp(-4*a));
I=sqrt(2*pi)*I;
end
