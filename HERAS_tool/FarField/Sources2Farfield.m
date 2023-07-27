classdef (Abstract) Sources2Farfield<POoperation
% Author : Francesco Lisi
% Revisions : v0.1.0 - 10/02/2023
% -------------------------------------------------------------------------
% Sources2Farfield is the parent class of the class Currents2Farfield
% that allows to compute the farfield  radiation patterns from the 
% currents on the reflector. 

    properties
        sources;
        coverage;
        D_RHCP_dB, D_LHCP_dB, D_V_dB, D_H_dB;
        E_Theta, E_Phi, E_H, E_V, E_RHCP, E_LHCP;
        copolar_component;
    end

    methods (Abstract)
        out=calculateFarField(obj,in);
    end

    methods (Access = protected)
      function cp = copyElement(obj)
         % Shallow copy object
         cp = copyElement@matlab.mixin.Copyable(obj);
         % Deep copy for the surfCurrents
         obj.sources = obj.sources.copy();
      end
    end

    methods 
        function [E_copolar,E_crosspolar,D_copolar_dB,D_crosspolar_dB]=ExtractCoPolXPolComponents(obj)
            switch obj.copolar_component
                case 'Theta'
                    E_copolar       = obj.E_Theta;
                    E_crosspolar    = obj.E_Phi;  
                    D_copolar_dB    = obj.D_Theta_dB;
                    D_crosspolar_dB = obj.D_Phi_dB;  
                case 'Phi'
                    E_copolar       = obj.E_Phi;
                    E_crosspolar    = obj.E_Theta;  
                    D_copolar_dB    = obj.D_Phi_dB;
                    D_crosspolar_dB = obj.D_Theta_dB;  
                case 'V'
                    E_copolar       = obj.E_V;
                    E_crosspolar    = obj.E_H;  
                    D_copolar_dB    = obj.D_V_dB;
                    D_crosspolar_dB = obj.D_H_dB;  
                case 'H'
                    E_copolar       = obj.E_H;
                    E_crosspolar    = obj.E_V;  
                    D_copolar_dB    = obj.D_H_dB;
                    D_crosspolar_dB = obj.D_V_dB;  
                case 'RHCP'
                    E_copolar       = obj.E_RHCP;
                    E_crosspolar    = obj.E_LHCP;  
                    D_copolar_dB    = obj.D_RHCP_dB;
                    D_crosspolar_dB = obj.D_LHCP_dB;  
                case 'LHCP'
                    E_copolar       = obj.E_LHCP;
                    E_crosspolar    = obj.E_RHCP;  
                    D_copolar_dB    = obj.D_LHCP_dB;
                    D_crosspolar_dB = obj.D_RHCP_dB;  
            end
        end

        function ApplyDistance(obj)

            distance_factor=exp(-1i*obj.sources.k_0*obj.coverage.r)./obj.coverage.r;

            % apply distance_factor
            obj.E_Theta = obj.E_Theta   .*distance_factor;
            obj.E_Phi   = obj.E_Phi     .*distance_factor;
            obj.E_H     = obj.E_H       .*distance_factor;
            obj.E_V     = obj.E_V       .*distance_factor;
            obj.E_RHCP  = obj.E_RHCP    .*distance_factor;
            obj.E_LHCP  = obj.E_LHCP    .*distance_factor;
            

        end

    end
end