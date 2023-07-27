%% Input
% uTheta, vPhi:             Far-field coordinates (either spherical theta/phi or sine-space uv)
% rx, phiy, z:              Reflector coordinates in reflector CS
% N:                        Jacobian
% J_x, J_y, J_z:            Equivalent electric surface currents in reflector CS
% k_0:                      Wavenumber
% rho_s:                    R-coordinate of reflector points in feed CS, only needed for Ludwig method
% phase_yn:                 Choice phase (depending on trapz or Ludwig method)
% uv_ThetaPhi:              Choice uv-coordinates or theta/phi-coordinates
% trapzLudwig:              Choice trapz or Ludwig method
% choice_GR:                Choice reflector grid
% d:                        Reflector diameter
% rx_single, phiy_single:   discrete x/y or rho/phi coordinates of reflector in reflector CS

%% Output
% E_x, E_y, E_z             Scattered electric field
% H_x, H_y, H_z             Scattered magnetic field            

%% Notes
% - With the help of the surface Jacobian transformation, the 3D integral
%   over the reflector surface can be reduced to an integration over the
%   projected circular region. Although integration is performed on the
%   projected aperture, the current is evaluated on the reflector surface
%   and the method is therefore different from the aperture distribution method
function [E_x, E_y, E_z, H_x, H_y, H_z] = NFIntegral(surfCurrents,cut)   
    
    %% Transform variables
    reflector = surfCurrents.reflector;
    X_FR=cut.X;
    Y_FR=cut.Y;
    Z_FR=cut.Z;
    [X_r, Y_r, Z_r] = reflector.getMeshPointsI('local');
    [rxG,phiyG,~] = reflector.getGridPoints();
    xrho_r = rxG(1,:);
    yphi_r = phiyG(:,1);
    [~,~,~,N] = reflector.getMeshNormals('local');
    J_x=surfCurrents.J_x;
    J_y=surfCurrents.J_y;
    J_z=surfCurrents.J_z;
    k0=surfCurrents.k_0;
    eta0=surfCurrents.eta_0;
    
%     eta0 = 120*pi;
    
    J_x = J_x .* -N;  % Multiply currents by Jacobian and norm. radii for 2D integration
    J_y = J_y .* -N;
    J_z = J_z .* -N;  
    
    [N_y,N_x,N_z] = size(X_FR);
    Int_Ex = zeros(N_y,N_x,N_z);
    Int_Ey = zeros(N_y,N_x,N_z);
    Int_Ez = zeros(N_y,N_x,N_z);
    Int_Hx = zeros(N_y,N_x,N_z);
    Int_Hy = zeros(N_y,N_x,N_z);
    Int_Hz = zeros(N_y,N_x,N_z);
    parfor l=1:N_y*N_x*N_z
        X = (X_FR(l) - X_r);
        Y = (Y_FR(l) - Y_r);
        Z = (Z_FR(l) - Z_r);
        [Phi,Theta,R]=cart2sph(X,Y,Z);

        %% Integration for E field
        G1 = (-1 - 1i*k0*R + k0^2*R.^2)./(R.^3);
        G2 = (3 + 3i*k0*R - k0^2*R.^2)./(R.^5);

        funx = (G1 .* J_x + (X_FR(l)-X_r).*G2 .* ((X_FR(l)-X_r).*J_x + (Y_FR(l)-Y_r).*J_y + (Z_FR(l)-Z_r).*J_z)) .* exp(-1i*k0*R);
        funy = (G1 .* J_y + (Y_FR(l)-Y_r).*G2 .* ((X_FR(l)-X_r).*J_x + (Y_FR(l)-Y_r).*J_y + (Z_FR(l)-Z_r).*J_z)) .* exp(-1i*k0*R);
        funz = (G1 .* J_z + (Z_FR(l)-Z_r).*G2 .* ((X_FR(l)-X_r).*J_x + (Y_FR(l)-Y_r).*J_y + (Z_FR(l)-Z_r).*J_z)) .* exp(-1i*k0*R);

        Int_Ex(l) = trapz(yphi_r,trapz(xrho_r,funx,2));
        Int_Ey(l) = trapz(yphi_r,trapz(xrho_r,funy,2));
        Int_Ez(l) = trapz(yphi_r,trapz(xrho_r,funz,2));  

        %% Integration for H field
        funx = (abs(Z_FR(l)-Z_r).*J_y - abs(Y_FR(l)-Y_r).*J_z) .* (1+1i*k0*R)./(R.^3) .* exp(-1i*k0*R);
        funy = (abs(X_FR(l)-X_r).*J_z - abs(Z_FR(l)-Z_r).*J_x) .* (1+1i*k0*R)./(R.^3) .* exp(-1i*k0*R);
        funz = (abs(Y_FR(l)-Y_r).*J_x - abs(X_FR(l)-X_r).*J_y) .* (1+1i*k0*R)./(R.^3) .* exp(-1i*k0*R);

        Int_Hx(l) = trapz(yphi_r,trapz(xrho_r,funx,2));
        Int_Hy(l) = trapz(yphi_r,trapz(xrho_r,funy,2));
        Int_Hz(l) = trapz(yphi_r,trapz(xrho_r,funz,2)); 
    end
    E_x = -1i*eta0/(4*pi*k0) * Int_Ex; 
    E_y = -1i*eta0/(4*pi*k0) * Int_Ey; 
    E_z = -1i*eta0/(4*pi*k0) * Int_Ez;
    
    H_x = 1/(4*pi) * Int_Hx; 
    H_y = 1/(4*pi) * Int_Hy; 
    H_z = 1/(4*pi) * Int_Hz;
end
