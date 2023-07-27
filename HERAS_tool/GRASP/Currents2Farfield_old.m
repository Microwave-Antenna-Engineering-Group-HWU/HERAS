classdef Currents2Farfield_old<Sources2Farfield_old
    % Currents2Farfield defines the cuts where the far field will be
    % calculated based on a surface current distribution.
    % First, create the object with Currents2Farfield(surfCurrents). Then,
    % define the cuts either with defineCuts (legacy code from Salva), or
    % directly dive the u,v points where the farfield values must be
    % calculated with definePlaneCuts (vectorized function). Finally, use
    % calculateFarField() to calculate the far field : it will call the
    % reflector associated with the surface currents object, reflector
    % which itself will call its mesh to do the integration.
    % Authors : Salvador Mercader Pellicer, Louis Dufour
    % Revisions : v0.1.0 - 25/10/2018

    properties
        surfCurrents;
    end
    methods (Access = protected)
      function cp = copyElement(obj)
         % Shallow copy object
         cp = copyElement@matlab.mixin.Copyable(obj);
         % Deep copy for the surfCurrents
         obj.surfCurrents = obj.surfCurrents.copy();
      end
   end
    methods
        function obj = Currents2Farfield_old()
            
        end
                    
        
        function calculateFarField(obj,surfCurrents,cut,trapzLudwig)
            obj.surfCurrents=surfCurrents;
            obj.cut=cut;
            if isa(obj.cut,'FarfieldSphericalGrid')
                uv_cut=obj.cut;
            elseif isa(obj.cut,'FarfieldSphericalCut')
                uv_cut=obj.cut.UV_Grid;
            end

            %% Farfield computation on CPU
            [obj.D_V_dB, obj.D_H_dB,obj.D_RHCP_dB, obj.D_LHCP_dB,...
                obj.E_Theta, obj.E_Phi, obj.E_V, obj.E_H, obj.E_RHCP, obj.E_LHCP] = ...
                integrateFFutility(obj.surfCurrents, uv_cut,trapzLudwig);


        end
    end
end

%% Auxiliary functions

function [D_V_dB, D_H_dB,D_RHCP_dB, D_LHCP_dB,...
    E_Theta, E_Phi, E_V, E_H, E_RHCP, E_LHCP] =  integrateFFutility(currents, cut,trapzLudwig)

reflector = currents.reflector;
type = reflector.mesh.intType;  % mesh type 'cartesian','cylindrical', etc..
Prad = currents.Prad;           % Prad of the feed

J_x = currents.J_x(:);
J_y = currents.J_y(:);
J_z = currents.J_z(:);
k_0 = currents.k_0;
eta_0 = currents.eta_0;
phase_yn = currents.phase_yn(:);   % 1 if phase added when creating feed; 0 if phase added before performing the integral (0 always when using Ludwig Method)
rho_s = currents.rho_s(:);         % Distance of each point in the integration grid with respect to the feed position

u = cut.U;                      % u grid
v = cut.V;                      % v grid
w = cut.W;                      % w=sqrt(1-(u.^2+v.^2))=cos(theta)
sinPhi_uv = cut.sinPhi_uv;      % sin(phi)
sinTheta_uv = cut.sinTheta_uv;  % sin(theta)
cosPhi_uv = cut.cosPhi_uv;      % cos(phi)
cosTheta_uv=cut.cosTheta_uv;    % cos(theta)
[xI, yI, zI] = reflector.getMeshPointsI('local');
xI=xI(:);
yI=yI(:);
zI=zI(:);
[~,~,~,N] = reflector.getMeshNormals('local');      % N=sqrt(nx^2+ny^2+nz^2)=sqrt(1+d_x(f)^2+d_y(f)^2) is the norm of the normal vector to the surface z=f(x,y)
N=N(:);

% Extract integration parameters: weights and Jacobian
switch type
    case {'cartesian','cylindrical'}
        % rxG,phiyG, rx_single and phiy_single are used in the Ludwig
        % integration. weights and Jacobian are used in the trapz
        % integration.
        weights = reflector.mesh.weights(:);
        [rxG,phiyG,Jacobian] = reflector.getGridPoints();   % Get the x,y or r,phi points on the integration mesh and Jacobian
        Jacobian=Jacobian(:);
        rx_single = rxG(1,:);
        phiy_single = phiyG(:,1);
    case 'triGauss'
        weights = reflector.mesh.weights(:);
        Jacobian = 1;
    case 'PolarGRASP'
        weights = reflector.mesh.weights(:);
        [~,~,Jacobian] = reflector.getGridPoints();
        Jacobian=Jacobian(:);
end

% N_uTheta=length(u(1,:));
% N_vPhi=length(u(:,1));

[N_v,N_u]=size(u);
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

switch trapzLudwig
    case 'trapz'
        % N is the norm of the n vector obtained as n=cross(v1,v2) where v1
        % and v2 are the derivatives of a surface along the two coordinates
        % defining it (see surface integral chapter in a textbook)
        J_x_N=J_x.*N.*exp(-1i*k_0*(rho_s)*(1-phase_yn));
        J_y_N=J_y.*N.*exp(-1i*k_0*(rho_s)*(1-phase_yn));
        J_z_N=J_z.*N.*exp(-1i*k_0*(rho_s)*(1-phase_yn));
        switch type
%             case {'cartesian','cylindrical'}
%                 parfor n=1:N_vPhi*N_uTheta
%                     Q=exp(1i*k_0*(u(n)*rx+v(n)*phiy+w(n)*z));
%                     funx=J_x_N.*Q;
%                     funy=J_y_N.*Q;
%                     funz=J_z_N.*Q;
%                     FTx(n)=trapz(phiy_single,trapz(rx_single,funx.*Jacobian,2));
%                     FTy(n)=trapz(phiy_single,trapz(rx_single,funy.*Jacobian,2));
%                     FTz(n)=trapz(phiy_single,trapz(rx_single,funz.*Jacobian,2));
%                 end
            case {'cartesian','cylindrical','triGauss','PolarGRASP'}
                parfor n=1:N_not_NaN
                    Q=exp(1i*k_0*(u(n)*xI+v(n)*yI+w(n)*zI));
                    funx=J_x_N.*Q;
                    funy=J_y_N.*Q;
                    funz=J_z_N.*Q;

                    %FTx(n,i) = sum(sum(obj.weights.*funx));
                    %FTy(n,i) = sum(sum(obj.weights.*funy));
                    %FTz(n,i) = sum(sum(obj.weights.*funz));
                    % slightly faster :
                    FTx_aux(n) = weights.'*(funx(:).*Jacobian);
                    FTy_aux(n) = weights.'*(funy(:).*Jacobian);
                    FTz_aux(n) = weights.'*(funz(:).*Jacobian);
                end
        end
    case 'ludwig'
        switch type
            case 'cartesian'
                [rx,phiy]=meshgrid(rx_single,phiy_single);
                N_y=length(J_x(:,1));
                N_x=length(J_x(1,:));
                delta_x=rx(1:N_y-1,2:N_x)-rx(1:N_y-1,1:N_x-1);
                delta_y=phiy(2:N_y,1:N_x-1)-phiy(1:N_y-1,1:N_x-1);
                F_x=J_x.*N;
                F_y=J_y.*N;
                F_z=J_z.*N;
                a_m_n_x=1/4*(3*F_x(1:N_y-1,1:N_x-1)-F_x(2:N_y,2:N_x)+F_x(2:N_y,1:N_x-1)+F_x(1:N_y-1,2:N_x));
                a_m_n_y=1/4*(3*F_y(1:N_y-1,1:N_x-1)-F_y(2:N_y,2:N_x)+F_y(2:N_y,1:N_x-1)+F_y(1:N_y-1,2:N_x));
                a_m_n_z=1/4*(3*F_z(1:N_y-1,1:N_x-1)-F_z(2:N_y,2:N_x)+F_z(2:N_y,1:N_x-1)+F_z(1:N_y-1,2:N_x));
                b_m_n_x=1./(2*delta_x).*(F_x(1:N_y-1,2:N_x)-F_x(1:N_y-1,1:N_x-1)+F_x(2:N_y,2:N_x)-F_x(2:N_y,1:N_x-1));
                b_m_n_y=1./(2*delta_x).*(F_y(1:N_y-1,2:N_x)-F_y(1:N_y-1,1:N_x-1)+F_y(2:N_y,2:N_x)-F_y(2:N_y,1:N_x-1));
                b_m_n_z=1./(2*delta_x).*(F_z(1:N_y-1,2:N_x)-F_z(1:N_y-1,1:N_x-1)+F_z(2:N_y,2:N_x)-F_z(2:N_y,1:N_x-1));
                c_m_n_x=1./(2*delta_y).*(F_x(2:N_y,1:N_x-1)-F_x(1:N_y-1,1:N_x-1)+F_x(2:N_y,2:N_x)-F_x(1:N_y-1,2:N_x));
                c_m_n_y=1./(2*delta_y).*(F_y(2:N_y,1:N_x-1)-F_y(1:N_y-1,1:N_x-1)+F_y(2:N_y,2:N_x)-F_y(1:N_y-1,2:N_x));
                c_m_n_z=1./(2*delta_y).*(F_z(2:N_y,1:N_x-1)-F_z(1:N_y-1,1:N_x-1)+F_z(2:N_y,2:N_x)-F_z(1:N_y-1,2:N_x));

                delta_x=delta_x(1,1); % This is true for equidistant grids
                delta_y=delta_y(1,1);
                jk=1i*k_0;
                lim=0.0001/k_0;
                parfor n=1:N_not_NaN
                    gamma=u(n)*rx+v(n)*phiy+z*w(n)-rho_s;
                    alpha_m_n=1/4*(3*gamma(1:N_y-1,1:N_x-1)-gamma(2:N_y,2:N_x)+gamma(2:N_y,1:N_x-1)+gamma(1:N_y-1,2:N_x));
                    beta_m_n=1./(2*delta_x).*(gamma(1:N_y-1,2:N_x)-gamma(1:N_y-1,1:N_x-1)+gamma(2:N_y,2:N_x)-gamma(2:N_y,1:N_x-1));
                    epsilon_m_n=1./(2*delta_y).*(gamma(2:N_y,1:N_x-1)-gamma(1:N_y-1,1:N_x-1)+gamma(2:N_y,2:N_x)-gamma(1:N_y-1,2:N_x));
                    expalpha=exp(jk*alpha_m_n);
                    expbeta=exp(jk*beta_m_n*delta_x);
                    jkbeta=jk*beta_m_n;
                    expepsilon=exp(jk*epsilon_m_n*delta_y);
                    jkepsilon=jk*epsilon_m_n;
                    aux11=(expbeta-1)./jkbeta;
                    aux12=(expepsilon-1)./jkepsilon;
                    aux1=aux11.*aux12;
                    aux21=(delta_x./jkbeta.*expbeta-(expbeta-1)./jkbeta.^2);
                    aux2=aux21.*aux12;
                    aux32=delta_y./jkepsilon.*expepsilon-(expepsilon-1)./jkepsilon.^2;
                    aux3=aux11.*aux32;
                    first=a_m_n_x.*aux1;
                    second=b_m_n_x.*aux2;
                    third=c_m_n_x.*aux3;
                    I_m_n_x=expalpha.*(first+second+third);
                    first=a_m_n_y.*aux1;
                    second=b_m_n_y.*aux2;
                    third=c_m_n_y.*aux3;
                    I_m_n_y=expalpha.*(first+second+third);
                    first=a_m_n_z.*aux1;
                    second=b_m_n_z.*aux2;
                    third=c_m_n_z.*aux3;
                    I_m_n_z=expalpha.*(first+second+third);

                    cond1=abs(epsilon_m_n(:));
                    cond2=abs(beta_m_n(:));
                    I_m_n_x(cond1<=lim & cond2>lim)=delta_y*expalpha(cond1<=lim & cond2>lim).*(aux11(cond1<=lim & cond2>lim).*(a_m_n_x(cond1<=lim & cond2>lim)+delta_y/2*c_m_n_x(cond1<=lim & cond2>lim))+b_m_n_x(cond1<=lim & cond2>lim).*aux21(cond1<=lim & cond2>lim));
                    I_m_n_y(cond1<=lim & cond2>lim)=delta_y*expalpha(cond1<=lim & cond2>lim).*(aux11(cond1<=lim & cond2>lim).*(a_m_n_y(cond1<=lim & cond2>lim)+delta_y/2*c_m_n_y(cond1<=lim & cond2>lim))+b_m_n_y(cond1<=lim & cond2>lim).*aux21(cond1<=lim & cond2>lim));
                    I_m_n_z(cond1<=lim & cond2>lim)=delta_y*expalpha(cond1<=lim & cond2>lim).*(aux11(cond1<=lim & cond2>lim).*(a_m_n_z(cond1<=lim & cond2>lim)+delta_y/2*c_m_n_z(cond1<=lim & cond2>lim))+b_m_n_z(cond1<=lim & cond2>lim).*aux21(cond1<=lim & cond2>lim));
                    I_m_n_x(cond1>lim & cond2<=lim)=delta_x*expalpha(cond1>lim & cond2<=lim).*(aux12(cond1>lim & cond2<=lim).*(a_m_n_x(cond1>lim & cond2<=lim)+delta_x/2*b_m_n_x(cond1>lim & cond2<=lim))+c_m_n_x(cond1>lim & cond2<=lim).*aux32(cond1>lim & cond2<=lim));
                    I_m_n_y(cond1>lim & cond2<=lim)=delta_x*expalpha(cond1>lim & cond2<=lim).*(aux12(cond1>lim & cond2<=lim).*(a_m_n_y(cond1>lim & cond2<=lim)+delta_x/2*b_m_n_y(cond1>lim & cond2<=lim))+c_m_n_y(cond1>lim & cond2<=lim).*aux32(cond1>lim & cond2<=lim));
                    I_m_n_z(cond1>lim & cond2<=lim)=delta_x*expalpha(cond1>lim & cond2<=lim).*(aux12(cond1>lim & cond2<=lim).*(a_m_n_z(cond1>lim & cond2<=lim)+delta_x/2*b_m_n_z(cond1>lim & cond2<=lim))+c_m_n_z(cond1>lim & cond2<=lim).*aux32(cond1>lim & cond2<=lim));
                    I_m_n_x(cond1<=lim & cond2<=lim)=delta_x*delta_y*expalpha(cond1<=lim & cond2<=lim).*(a_m_n_x(cond1<=lim & cond2<=lim)+delta_x/2*b_m_n_x(cond1<=lim & cond2<=lim)+delta_y/2*c_m_n_x(cond1<=lim & cond2<=lim));
                    I_m_n_y(cond1<=lim & cond2<=lim)=delta_x*delta_y*expalpha(cond1<=lim & cond2<=lim).*(a_m_n_y(cond1<=lim & cond2<=lim)+delta_x/2*b_m_n_y(cond1<=lim & cond2<=lim)+delta_y/2*c_m_n_y(cond1<=lim & cond2<=lim));
                    I_m_n_z(cond1<=lim & cond2<=lim)=delta_x*delta_y*expalpha(cond1<=lim & cond2<=lim).*(a_m_n_z(cond1<=lim & cond2<=lim)+delta_x/2*b_m_n_z(cond1<=lim & cond2<=lim)+delta_y/2*c_m_n_z(cond1<=lim & cond2<=lim));

                    FTx_aux(n)=sum(sum(I_m_n_x,2));
                    FTy_aux(n)=sum(sum(I_m_n_y,2));
                    FTz_aux(n)=sum(sum(I_m_n_z,2));
                end
            case 'cylindrical'
                [rx,phiy]=meshgrid(rx_single,phiy_single);
                F_x=J_x.*N.*rx;
                F_y=J_y.*N.*rx;
                F_z=J_z.*N.*rx;
                N_phi=length(J_x(:,1));
                N_r=length(J_x(1,:));
                a_m_n_x=1/4*(3*F_x(1:N_phi-1,1:N_r-1)-F_x(2:N_phi,2:N_r)+F_x(2:N_phi,1:N_r-1)+F_x(1:N_phi-1,2:N_r));
                a_m_n_y=1/4*(3*F_y(1:N_phi-1,1:N_r-1)-F_y(2:N_phi,2:N_r)+F_y(2:N_phi,1:N_r-1)+F_y(1:N_phi-1,2:N_r));
                a_m_n_z=1/4*(3*F_z(1:N_phi-1,1:N_r-1)-F_z(2:N_phi,2:N_r)+F_z(2:N_phi,1:N_r-1)+F_z(1:N_phi-1,2:N_r));
                delta_r=rx(1:N_phi-1,2:N_r)-rx(1:N_phi-1,1:N_r-1);
                b_m_n_x=1./(2*delta_r).*(F_x(1:N_phi-1,2:N_r)-F_x(1:N_phi-1,1:N_r-1)+F_x(2:N_phi,2:N_r)-F_x(2:N_phi,1:N_r-1));
                b_m_n_y=1./(2*delta_r).*(F_y(1:N_phi-1,2:N_r)-F_y(1:N_phi-1,1:N_r-1)+F_y(2:N_phi,2:N_r)-F_y(2:N_phi,1:N_r-1));
                b_m_n_z=1./(2*delta_r).*(F_z(1:N_phi-1,2:N_r)-F_z(1:N_phi-1,1:N_r-1)+F_z(2:N_phi,2:N_r)-F_z(2:N_phi,1:N_r-1));
                delta_phi=phiy(2:N_phi,1:N_r-1)-phiy(1:N_phi-1,1:N_r-1);
                c_m_n_x=1./(2*delta_phi).*(F_x(2:N_phi,1:N_r-1)-F_x(1:N_phi-1,1:N_r-1)+F_x(2:N_phi,2:N_r)-F_x(1:N_phi-1,2:N_r));
                c_m_n_y=1./(2*delta_phi).*(F_y(2:N_phi,1:N_r-1)-F_y(1:N_phi-1,1:N_r-1)+F_y(2:N_phi,2:N_r)-F_y(1:N_phi-1,2:N_r));
                c_m_n_z=1./(2*delta_phi).*(F_z(2:N_phi,1:N_r-1)-F_z(1:N_phi-1,1:N_r-1)+F_z(2:N_phi,2:N_r)-F_z(1:N_phi-1,2:N_r));

                delta_r=delta_r(1,1); % This is true for equidistant grids
                delta_phi=delta_phi(1,1);
                cosphi=cos(phiy);
                sinphi=sin(phiy);
                jk=1i*k_0;
                lim=0.00001/k_0;
                parfor n=1:N_not_NaN
                    gamma=rx.*(u(n)*cosphi+v(n)*sinphi)+z*w(n)-rho_s;
                    alpha_m_n=1/4*(3*gamma(1:N_phi-1,1:N_r-1)-gamma(2:N_phi,2:N_r)+gamma(2:N_phi,1:N_r-1)+gamma(1:N_phi-1,2:N_r));
                    beta_m_n=1./(2*delta_r).*(gamma(1:N_phi-1,2:N_r)-gamma(1:N_phi-1,1:N_r-1)+gamma(2:N_phi,2:N_r)-gamma(2:N_phi,1:N_r-1));
                    epsilon_m_n=1./(2*delta_phi).*(gamma(2:N_phi,1:N_r-1)-gamma(1:N_phi-1,1:N_r-1)+gamma(2:N_phi,2:N_r)-gamma(1:N_phi-1,2:N_r));
                    expalpha=exp(jk*alpha_m_n);
                    expbeta=exp(jk*beta_m_n*delta_r);
                    jkbeta=jk*beta_m_n;
                    expepsilon=exp(jk*epsilon_m_n*delta_phi);
                    jkepsilon=jk*epsilon_m_n;
                    aux11=(expbeta-1)./jkbeta;
                    aux12=(expepsilon-1)./jkepsilon;
                    aux1=aux11.*aux12;
                    aux21=(delta_r./jkbeta.*expbeta-(expbeta-1)./jkbeta.^2);
                    aux2=aux21.*aux12;
                    aux32=delta_phi./jkepsilon.*expepsilon-(expepsilon-1)./jkepsilon.^2;
                    aux3=aux11.*aux32;
                    first=a_m_n_x.*aux1;
                    second=b_m_n_x.*aux2;
                    third=c_m_n_x.*aux3;
                    I_m_n_x=expalpha.*(first+second+third);
                    first=a_m_n_y.*aux1;
                    second=b_m_n_y.*aux2;
                    third=c_m_n_y.*aux3;
                    I_m_n_y=expalpha.*(first+second+third);
                    first=a_m_n_z.*aux1;
                    second=b_m_n_z.*aux2;
                    third=c_m_n_z.*aux3;
                    I_m_n_z=expalpha.*(first+second+third);

                    cond1=abs(epsilon_m_n(:));
                    cond2=abs(beta_m_n(:));
                    I_m_n_x(cond1<=lim & cond2>lim)=delta_phi*expalpha(cond1<=lim & cond2>lim).*(aux11(cond1<=lim & cond2>lim).*(a_m_n_x(cond1<=lim & cond2>lim)+delta_phi/2*c_m_n_x(cond1<=lim & cond2>lim))+b_m_n_x(cond1<=lim & cond2>lim).*aux21(cond1<=lim & cond2>lim));
                    I_m_n_y(cond1<=lim & cond2>lim)=delta_phi*expalpha(cond1<=lim & cond2>lim).*(aux11(cond1<=lim & cond2>lim).*(a_m_n_y(cond1<=lim & cond2>lim)+delta_phi/2*c_m_n_y(cond1<=lim & cond2>lim))+b_m_n_y(cond1<=lim & cond2>lim).*aux21(cond1<=lim & cond2>lim));
                    I_m_n_z(cond1<=lim & cond2>lim)=delta_phi*expalpha(cond1<=lim & cond2>lim).*(aux11(cond1<=lim & cond2>lim).*(a_m_n_z(cond1<=lim & cond2>lim)+delta_phi/2*c_m_n_z(cond1<=lim & cond2>lim))+b_m_n_z(cond1<=lim & cond2>lim).*aux21(cond1<=lim & cond2>lim));
                    I_m_n_x(cond1>lim & cond2<=lim)=delta_r*expalpha(cond1>lim & cond2<=lim).*(aux12(cond1>lim & cond2<=lim).*(a_m_n_x(cond1>lim & cond2<=lim)+delta_r/2*b_m_n_x(cond1>lim & cond2<=lim))+c_m_n_x(cond1>lim & cond2<=lim).*aux32(cond1>lim & cond2<=lim));
                    I_m_n_y(cond1>lim & cond2<=lim)=delta_r*expalpha(cond1>lim & cond2<=lim).*(aux12(cond1>lim & cond2<=lim).*(a_m_n_y(cond1>lim & cond2<=lim)+delta_r/2*b_m_n_y(cond1>lim & cond2<=lim))+c_m_n_y(cond1>lim & cond2<=lim).*aux32(cond1>lim & cond2<=lim));
                    I_m_n_z(cond1>lim & cond2<=lim)=delta_r*expalpha(cond1>lim & cond2<=lim).*(aux12(cond1>lim & cond2<=lim).*(a_m_n_z(cond1>lim & cond2<=lim)+delta_r/2*b_m_n_z(cond1>lim & cond2<=lim))+c_m_n_z(cond1>lim & cond2<=lim).*aux32(cond1>lim & cond2<=lim));
                    I_m_n_x(cond1<=lim & cond2<=lim)=delta_r*delta_phi*expalpha(cond1<=lim & cond2<=lim).*(a_m_n_x(cond1<=lim & cond2<=lim)+delta_r/2*b_m_n_x(cond1<=lim & cond2<=lim)+delta_phi/2*c_m_n_x(cond1<=lim & cond2<=lim));
                    I_m_n_y(cond1<=lim & cond2<=lim)=delta_r*delta_phi*expalpha(cond1<=lim & cond2<=lim).*(a_m_n_y(cond1<=lim & cond2<=lim)+delta_r/2*b_m_n_y(cond1<=lim & cond2<=lim)+delta_phi/2*c_m_n_y(cond1<=lim & cond2<=lim));
                    I_m_n_z(cond1<=lim & cond2<=lim)=delta_r*delta_phi*expalpha(cond1<=lim & cond2<=lim).*(a_m_n_z(cond1<=lim & cond2<=lim)+delta_r/2*b_m_n_z(cond1<=lim & cond2<=lim)+delta_phi/2*c_m_n_z(cond1<=lim & cond2<=lim));

                    FTx_aux(n)=sum(sum(I_m_n_x,2));
                    FTy_aux(n)=sum(sum(I_m_n_y,2));
                    FTz_aux(n)=sum(sum(I_m_n_z,2));
                end
            case {'triGauss','PolarGRASP'}
                error('Ludwig method not developed for triangular/PolarGRASP mesh');
        end
end

FTx(not_NaN_mask)=FTx_aux;
FTy(not_NaN_mask)=FTy_aux;
FTz(not_NaN_mask)=FTz_aux;
clear FTx_aux FTy_aux FTz_aux;

% w is cosTheta_uv
E_Theta=-1i*k_0*eta_0/(4*pi).*(cosTheta_uv.*cosPhi_uv.*FTx+cosTheta_uv.*sinPhi_uv.*FTy-sinTheta_uv.*FTz);
E_Phi=-1i*k_0*eta_0/(4*pi).*(-sinPhi_uv.*FTx+cosPhi_uv.*FTy);
E_Theta_norm=E_Theta;%/max(sqrt(abs(E_Theta_PO_1_1(:)).^2+abs(E_Phi_PO_1_1(:)).^2));
E_Phi_norm=E_Phi;%/max(sqrt(abs(E_Theta_PO_1_1(:)).^2+abs(E_Phi_PO_1_1(:)).^2));
E_V=E_Theta_norm.*cosPhi_uv-E_Phi_norm.*sinPhi_uv;
E_H=E_Theta_norm.*sinPhi_uv+E_Phi_norm.*cosPhi_uv;
D_V_dB=10*log10(abs(E_V).^2/(2*eta_0)*4*pi/Prad);
D_H_dB=10*log10(abs(E_H).^2/(2*eta_0)*4*pi/Prad);
E_RHCP=1/sqrt(2)*(E_V+1i*(E_H));
E_LHCP=1/sqrt(2)*(E_V-1i*(E_H));
D_RHCP_dB=10*log10(abs(E_RHCP).^2/(2*eta_0)*4*pi/Prad);
D_LHCP_dB=10*log10(abs(E_LHCP).^2/(2*eta_0)*4*pi/Prad); 
end