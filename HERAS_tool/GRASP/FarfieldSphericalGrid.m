classdef FarfieldSphericalGrid<POoperation & dynamicprops
    %FARFIELDSPHERICALGRID Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        min_u
        max_u
        N_u
        min_v
        max_v
        N_v
        u_single
        v_single
        U
        V
        W
        sinPhi_uv
        sinTheta_uv
        cosPhi_uv
        cosTheta_uv
    end
    
    methods
        function obj = FarfieldSphericalGrid()
            %FARFIELDSPHERICALCUT Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function defineCutsFromSingle(obj,min_u,max_u,N_u,min_v,max_v,N_v)
            %DEFINECUTS Summary of this method goes here
            %   Detailed explanation goes here
            obj.min_u=min_u;
            obj.max_u=max_u;
            obj.N_u=N_u;
            obj.min_v=min_v;
            obj.max_v=max_v;
            obj.N_v=N_v;
            obj.u_single=linspace(obj.min_u,obj.max_u,obj.N_u);
            obj.v_single=linspace(obj.min_v,obj.max_v,obj.N_v);
            [obj.U,obj.V]=meshgrid(obj.u_single,obj.v_single);
            obj.completeCuts();
        end
        
        function defineCutsFromMatrix(obj,U,V)
            %DEFINECUTS Summary of this method goes here
            %   Detailed explanation goes here
            obj.U=U;
            obj.V=V;
            obj.completeCuts();
        end
        function completeCuts(obj)
            %DEFINECUTS Summary of this method goes here
            %   Detailed explanation goes here
            obj.W=sqrt(1-((obj.U).^2+(obj.V).^2));
            obj.W(sqrt((obj.U).^2+(obj.V).^2)>1)=0;
            obj.sinPhi_uv=(obj.V)./sqrt((obj.U).^2+(obj.V).^2);
            obj.sinPhi_uv(sqrt((obj.U).^2+(obj.V).^2)>1)=0;
            obj.sinPhi_uv(sqrt((obj.U).^2+(obj.V).^2)==0)=0;
            obj.cosPhi_uv=(obj.U)./sqrt((obj.U).^2+(obj.V).^2);
            obj.cosPhi_uv(sqrt((obj.U).^2+(obj.V).^2)>1)=0;
            obj.cosPhi_uv(sqrt((obj.U).^2+(obj.V).^2)==0)=1;
            obj.sinTheta_uv=sqrt((obj.U).^2+(obj.V).^2);
            obj.cosTheta_uv=obj.W;
        end
        
    end
end

