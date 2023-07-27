classdef CylindricalMesh<POReflectorMesh
    % CylindricalMesh implements a regular cylindrical mesh on a reflector
    % with a disk aperture.
    % Integration is trapz, using Matlab bult in function and the parallel
    % processing toolbox
    % Authors : Salvador Mercader Pellicer, Louis Dufour
    % v0.1.0 : 25/10/2018
    properties
        Nrho, Nphi, lambda_0;
        rho_single, phi_single;
        rhoI, phiI;
        Jacobian;
        weights;
    end
    
    methods
        function obj = CylindricalMesh(Nrho, Nphi, lambda_0)
            obj.Nrho = Nrho;
            obj.Nphi = Nphi;
            obj.lambda_0 = lambda_0;
            obj.intType = 'cylindrical';
        end
        
        function meshNow(obj, rim, lines)
            switch rim.type
                case 'elliptic'
                    if rim.a1 ~= rim.a2
                        x0 = rim.x0;
                        y0 = rim.y0;
                        a1 = rim.a1;
                        a2 = rim.a2;
                        obj.rho_single=linspace(0,1,obj.Nrho);
                        obj.phi_single=linspace(0,2*pi,obj.Nphi);
                        [obj.rhoI,obj.phiI]=meshgrid(obj.rho_single,obj.phi_single);
                        obj.xI = x0 + a1 * obj.rhoI.*cos(obj.phiI);
                        obj.yI= y0 + a2 * obj.rhoI.*sin(obj.phiI);
                        obj.Jacobian = obj.rhoI*a1*a2;
                        obj.xE = obj.xI;
                        obj.yE = obj.yI;
                        w_rho_single=(obj.rho_single(end)-obj.rho_single(1))/(2*(obj.Nrho-1))*[1, 2*ones(1,obj.Nrho-2),1];
                        w_phi_single=(obj.phi_single(end)-obj.phi_single(1))/(2*(obj.Nphi-1))*[1, 2*ones(1,obj.Nphi-2),1];

                        [Wrho,Wphi]=meshgrid(w_rho_single,w_phi_single);
                        W=Wrho.*Wphi;
                        obj.weights=W;
                        obj.isMeshed = true;
                    else
                        d = 2*rim.a1;
                        x0 = rim.x0;
                        y0 = rim.y0;
                        obj.rho_single=linspace(0,d/2,obj.Nrho);
                        obj.phi_single=linspace(0,2*pi,obj.Nphi);
                        [obj.rhoI,obj.phiI]=meshgrid(obj.rho_single,obj.phi_single);
                        obj.xI = x0 + obj.rhoI.*cos(obj.phiI);
                        obj.yI= y0 + obj.rhoI.*sin(obj.phiI);
                        obj.Jacobian = obj.rhoI;
                        obj.xE = obj.xI;
                        obj.yE = obj.yI;
                        w_rho_single=(obj.rho_single(end)-obj.rho_single(1))/(2*(obj.Nrho-1))*[1, 2*ones(1,obj.Nrho-2),1];
                        w_phi_single=(obj.phi_single(end)-obj.phi_single(1))/(2*(obj.Nphi-1))*[1, 2*ones(1,obj.Nphi-2),1];

                        [Wrho,Wphi]=meshgrid(w_rho_single,w_phi_single);
                        W=Wrho.*Wphi;
                        obj.weights=W;
                        obj.isMeshed = true;
                    end
                otherwise
                    error('CartesianMesh only works with elliptic rims');
            end
        end
        
        function [x1,x2,J] = getGridPoints(obj)
            x1 = obj.rhoI;
            x2 = obj.phiI;
            J = obj.Jacobian;
        end
    end
end