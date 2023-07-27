classdef Coverage < POobject
% Author : Francesco Lisi
% Revisions : 0.1.0 - 06/03/2023
% -------------------------------------------------------------------------
% This class defines the coverage, i.e. the set of farfield points
% where the fields are computed. Each farfield point is defined by
% its theta and phi coordinates with respect to the coverage CS.
% Alternatively one can use the three coordinates of the associated unit
% vector, i.e. U=sin(Theta)cos(Phi), V=sin(Theta)sin(Phi), and
% W=cos(Theta).


    properties
        U, V, W;
        Theta, Phi;
        sinTheta, cosTheta, sinPhi, cosPhi;
        index = [];
        r;
    end

    methods
        function obj = Coverage(x,y,z,ex,ey,ez)
        % Constructor
            ex = ex/norm(ex); ey = ey/norm(ey); ez = ez/norm(ez);
            if dot(ex,ey)>1e-10 ||dot(ex,ez)>1e-10 ||dot(ey,ez)>1e-10
                warning('The coordinate syste defined for the feed is not orthogonal');
                return
            end
            obj.position = [x;y;z];
            obj.direction = [ex(:) ey(:) ez(:)];
        end

        function DefineCoverageFromUVW(obj,U,V,W)
            % This method defines the coverage from the u, v, and w values

            % Creation of a mask containing the points whose norm is equal
            % to 1 (within a tolerance defined by tolerance)
            tolerance=1e-9;
            norm=sqrt(U.^2+V.^2+W.^2);
            norm_equal_1_mask=abs(norm-1)<=tolerance;

            % Set to NaN points outside the mask
            U(~norm_equal_1_mask)=NaN;
            V(~norm_equal_1_mask)=NaN;
            W(~norm_equal_1_mask)=NaN;

            % Normalize the remaining points such that their norm is
            % exactly 1
            U=U./norm;
            V=V./norm;
            W=W./norm;

            % Store U, V, and W
            obj.U = U;
            obj.V = V;
            obj.W = W;

            % Compute Theta and Phi from U, V, and W
            [obj.Theta,obj.Phi]=UVWtoThetaPhi(U,V,W);

            % Set the remaining properties
            obj.CompleteCoverage();
        end

        function DefineCoverageFromUV(obj,U,V)
            % This method defines the coverage from the u, and v values.
            % The w value is assumed to be positive, so since w=cos(Theta),
            % Theta belongs to the interval [0,pi/2]

            % Compute mask of points outside unit circle
            mask=(U.^2+V.^2)>1;

            % Set to NaN points outside the unit circle
            U(mask)=NaN;
            V(mask)=NaN;

            % Compute W from U and V
            W=sqrt(1-(U.^2+V.^2));

            % Call DefineCoverageFromUVW method
            obj.DefineCoverageFromUVW(U,V,W);
        end

        function DefineCoverageFromThetaPhi(obj,Theta,Phi)
            % This method defines the coverage from the Theta and Phi
            % coordinates

            % Store Theta and Phi
            obj.Theta=Theta;
            obj.Phi=Phi;

            % Compute U, V, an W from Theta and Phi
            [obj.U,obj.V,obj.W]=ThetaPhitoUVW(Theta,Phi);

            % Set the remaining properties
            obj.CompleteCoverage();
        end

        function CompleteCoverage(obj)
            % This method computes the remaining properties

            % Compute sin and cosine of Theta and Phi
            obj.sinTheta=sin(obj.Theta);
            obj.cosTheta=cos(obj.Theta);
            obj.sinPhi=sin(obj.Phi);
            obj.cosPhi=cos(obj.Phi);

            % set index in case the index has not been set previously
            if isempty(obj.index)
                [N,M]=size(obj.Theta);
                [column_index,row_index]=meshgrid(1:M,1:N);
                obj.index = [row_index(:),column_index(:)];
            end
        end

        function ComputeDistance(obj,earthCS,earth_radius)
            % This method computes the distance of each coverage point from
            % the coverage CS

            % coverage point unit vector in the coverage CS
            rFF_unit=[obj.U(:).';obj.V(:).';obj.W(:).'];

            % b is the vector from the Earth CS origin to the coverage CS
            % origin expressed in the global CS
            b=obj.position-earthCS.position;

            % obj.direction.'*b corresponds to the vector b expressed in
            % the coverage CS. So v is the inner product between each FF
            % unit vector and b.
            v=rFF_unit.'*obj.direction.'*b;

            % r is the distance correponding to the intersection between
            % the Earth surface and the line passing through the ceverage
            % CS origin with direction rFF_unit
            r=-v-sqrt(v.^2+earth_radius^2-norm(b)^2);
            
            % assign distance to coverage
            obj.r=reshape(r,size(obj.U));

        end

        function disp2fig(obj,ax, CSref)
            % Display Coverage CS

            if nargin == 1
                ax = [];
            elseif nargin == 2
                CSref = CSlocal(0,0,0,[1,0,0],[0,1,0],[0,0,1]);
            end

            if isempty(ax)
                ax = figure;
            else
                if isa(ax,'matlab.graphics.axis.Axes')
                    axes(ax);  hold on;
                elseif isa(ax,'matlab.ui.control.UIAxes')
                    hold(ax,'on');
                else
                    figure(ax); hold on;
                end
            end


            % Triade
            obj.plotTriade(ax, CSref);
        end
        
    end
end