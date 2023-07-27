classdef PolarMeshGRASP<POReflectorMesh
    % PolarMeshGRASP implements an irregular polar mesh on a reflector
    % with a disk aperture.
    % Integration is trapz for phi component and the rule of Gauss for rho
    % component, using weights and the parallel processing toolbox
    % A detailed explanation of the irregular polar grid can be found in
    % the document 'GRASP_Technical_Description' paragraph 3.1.2.1. The
    % document is accessible by downloading the GRASP software by TICRA.
    % Authors   : Alejandro Baldominos
    % v0.1.0    : 01/11/2021
    properties
        Nrho, Nphi;
        rho_single, Nphi_single;
        rhoI, phiI;
        tri;
        Jacobian;
        weights;
    end
    
    methods
        function obj = PolarMeshGRASP(Nrho, Nphi)
            obj.Nrho    = Nrho;
            obj.Nphi    = Nphi;
            obj.intType = 'PolarGRASP';
        end
        
        function meshNow(obj, rim, lines)
            switch rim.type
                case 'elliptic'
                    if rim.a1 ~= rim.a2                        
                        x0 = rim.x0;
                        y0 = rim.y0;
                        a1 = rim.a1;
                        a2 = rim.a2;
                        
                        [Rho,Phi,W] = obj.PolarMesh(1);

                        [X,Y] = pol2cart(Phi,Rho);

                        obj.rhoI = Rho(:);
                        obj.phiI = Phi(:);

                        obj.xI          = x0 + a1 * X(:);
                        obj.yI          = y0 + a2 * Y(:);
                        obj.Jacobian    = obj.rhoI*a1*a2;
                        obj.weights     = W(:);
                        obj.xE          = obj.xI;
                        obj.yE          = obj.yI;
                        obj.tri         = delaunay(obj.xE,obj.yE);
                        obj.isMeshed    = true;
                    else
                        d   = 2*rim.a1;
                        x0  = rim.x0;
                        y0  = rim.y0;
                        
                        [Rho,Phi,W] = obj.PolarMesh(d);

                        [X,Y] = pol2cart(Phi,Rho);

                        
                        obj.rhoI = Rho(:);
                        obj.phiI = Phi(:);

                        obj.xI          = x0 + X(:);
                        obj.yI          = y0 + Y(:);
                        obj.Jacobian    = obj.rhoI;
                        obj.weights     = W(:);
                        obj.xE          = obj.xI;
                        obj.yE          = obj.yI;
                        obj.tri         = delaunay(obj.xE,obj.yE);
                        obj.isMeshed    = true;
                    end
                otherwise
                    error('PolarMeshGRASP only works with elliptic rims');
            end
        end
        
        function [Rho,Phi,W] = PolarMesh(obj,d)
            [obj.rho_single,W_rho_single]   = lgwt(obj.Nrho,0,d/2);
            obj.Nphi_single                 = round((max([obj.Nphi,10])-10)*sqrt(obj.rho_single./max(obj.rho_single)))+10;
            Rho                             = zeros(1,sum(obj.Nphi_single));
            W_Rho                           = Rho;
            Phi                             = Rho;
            W_Phi                           = Rho;

            j=1;
            for i=1:obj.Nrho

                Rho(j:j+obj.Nphi_single(i)-1)   = obj.rho_single(i)*ones(1,obj.Nphi_single(i));
                W_Rho(j:j+obj.Nphi_single(i)-1) = W_rho_single(i)*ones(1,obj.Nphi_single(i));
                Phi_aux                         = linspace(0,2*pi,obj.Nphi_single(i));
                Phi(j:j+obj.Nphi_single(i)-1)   = Phi_aux(1:end);
                W_Phi(j:j+obj.Nphi_single(i)-1) = (2*pi/(2*(obj.Nphi_single(i)-1)))*[1, 2*ones(1,obj.Nphi_single(i)-2),1];
                j                               = j + obj.Nphi_single(i);
            
            end
            W=W_Rho.*W_Phi;
        end
        
        function [x1,x2,J] = getGridPoints(obj)
            x1  = obj.rhoI;
            x2  = obj.phiI;
            J   = obj.Jacobian;
        end
    end
end