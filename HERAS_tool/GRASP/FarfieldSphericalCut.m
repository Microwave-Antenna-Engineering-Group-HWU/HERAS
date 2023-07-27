classdef FarfieldSphericalCut<POoperation
    %FARFIELDSPHERICALCUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        min_Theta
        max_Theta
        N_Theta
        min_Phi
        max_Phi
        N_Phi
        Theta_single
        Phi_single
        Theta
        Phi
        UV_Grid
    end
    
    methods
        function obj = FarfieldSphericalCut()
            %FARFIELDSPHERICALCUT Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function defineCutsFromSingle(obj,min_Theta,max_Theta,N_Theta,min_Phi,max_Phi,N_Phi)
            %DEFINECUTS Summary of this method goes here
            %   Detailed explanation goes here
            obj.min_Theta=min_Theta;
            obj.max_Theta=max_Theta;
            obj.N_Theta=N_Theta;
            obj.min_Phi=min_Phi;
            obj.max_Phi=max_Phi;
            obj.N_Phi=N_Phi;
            obj.Theta_single=linspace(obj.min_Theta,obj.max_Theta,obj.N_Theta);
            obj.Phi_single=linspace(obj.min_Phi,obj.max_Phi,obj.N_Phi);
            [obj.Theta,obj.Phi]=meshgrid(obj.Theta_single,obj.Phi_single);
            obj.completeCuts();
        end
        
        function defineCutsFromMatrix(obj,Theta,Phi)
            %DEFINECUTS Summary of this method goes here
            %   Detailed explanation goes here
            obj.Theta=Theta;
            obj.Phi=Phi;
            obj.completeCuts();
        end
        
         function completeCuts(obj)
            %DEFINECUTS Summary of this method goes here
            %   Detailed explanation goes here
            U=sin(obj.Theta).*cos(obj.Phi);
            V=sin(obj.Theta).*sin(obj.Phi);
            obj.UV_Grid=FarfieldSphericalGrid();
            obj.UV_Grid.defineCutsFromMatrix(U,V);
            obj.UV_Grid.sinPhi_uv=sin(obj.Phi);
            obj.UV_Grid.sinTheta_uv=sin(obj.Theta);
            obj.UV_Grid.cosPhi_uv=cos(obj.Phi);
            obj.UV_Grid.cosTheta_uv=cos(obj.Theta);
            obj.UV_Grid.W=cos(obj.Theta);
        end
        
    end
end

