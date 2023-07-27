classdef (Abstract) Sources2Farfield_old<POoperation
    % Sources2Farfield is the parent class of the classes Currents2Farfield
    % that allow to compute the farfield radiation patterns from the 
    % currents on the reflector 
    % Authors : Francesco Lisi
    % Revisions : v0.1.0 - 10/02/2023
    properties
        cut;
        %Phi_uv, Theta_uv, u, v, w, N_u, N_v, cosPhi_uv, sinPhi_uv, sinTheta_uv;
        D_RHCP_dB, D_LHCP_dB, D_V_dB, D_H_dB;
        E_Theta, E_Phi, E_H, E_V, E_RHCP, E_LHCP;
    end

    methods (Abstract)
        out=calculateFarField(obj,in);
    end
end