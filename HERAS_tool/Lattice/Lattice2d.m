classdef Lattice2d < Lattice
    % This class defines a general 2D lattice.
    %
    % Author : Francesco Lisi
    % Revisions : 0.1.0 - 28/02/2023

    properties
        rotation_angle = 0;
    end
 
    methods
        function obj=Lattice2d(unit_cell)
            % Constructor class. One colour is assumed by default.
            obj@Lattice(unit_cell);
            obj.N_colors=1;
            obj.M_colors_matrix=eye(obj.N_dimensions);
            obj.color_step_translation=zeros(obj.N_dimensions,obj.N_colors);
        end

        function Rotate(obj,rotation_angle)
            % This method allows to rotate the 2D lattice. The sense of
            % rotatio n is counter clockwise.
            obj.rotation_angle=rotation_angle;
            rotation_matrix=[cos(rotation_angle),-sin(rotation_angle);sin(rotation_angle),cos(rotation_angle)];
            obj.unit_cell = rotation_matrix*obj.unit_cell;
        end

        function SetMultipleColors(obj,N_colors)
            % Assign N_colors property
            obj.N_colors=N_colors;

            % Assign M_colors_matrix depending on the number of colors
            switch N_colors
                case 1
                     obj.M_colors_matrix=eye(2,2);                   
                case 3
                    obj.M_colors_matrix=[1, 2; 1 , -1];
                case 4
                    obj.M_colors_matrix=[2,0;0,2];
                case 7
                    obj.M_colors_matrix=[2 , 3 ; 1 , -2];
                otherwise
                    error('This functions supports only 3,4 or 7 colors.')
            end

            %% Computation of each colour translation vector
            % fundamental_parallelogram_minusDelta contains the four 
            % vertexes of the fundamental parallelogram slightly reduced to 
            % exclude points on the perimeter.
            % fundamental_parallelogram_half_perimeter contains the four
            % vertexes of a polygon that only includes half perimeter of
            % the fundamental parallelogram.
            delta=eps;
            fundamental_parallelogram_minusDelta=       obj.M_colors_matrix*[1-delta , 0+delta, 0+delta, 1-delta ; 0+delta , 0+delta, 1-delta , 1-delta];
            fundamental_parallelogram_half_perimeter=   obj.M_colors_matrix*[1-delta , 0-delta, 0, 0+delta ; 0 , 0-delta, 1-delta , 0+delta];
            
            % The following lines of code compute all the translation vectors defined
            % as the vectors n in N^2 that satisfy the following property: n belongs to
            % the fundamental parallelogram defined by M_colors_matrix, i.e. the set of
            % points r in R^2 expressed as r=M_colors_matrix*s for all s in [0,1)^2.
            % Check the paper "MULTIPLE BEAMS FROM PLANAR ARRAYS" by Piero Angeletti
            N_max=sum(abs(obj.M_colors_matrix*ones(2,1)));
            translation_vectors=[kron((-N_max:N_max),ones(1,2*N_max+1));kron(ones(1,2*N_max+1),(-N_max:N_max))];
            translation_vectors_inside_parallelogram_minusDelta_index=inpolygon(translation_vectors(1,:),translation_vectors(2,:),fundamental_parallelogram_minusDelta(1,:),fundamental_parallelogram_minusDelta(2,:));
            translation_vectors_inside_parallelogram_half_perimeter_index=inpolygon(translation_vectors(1,:),translation_vectors(2,:),fundamental_parallelogram_half_perimeter(1,:),fundamental_parallelogram_half_perimeter(2,:));
            translation_vectors=[translation_vectors(:,translation_vectors_inside_parallelogram_minusDelta_index),translation_vectors(:,translation_vectors_inside_parallelogram_half_perimeter_index)];
            [~,sort_index]=sort(vecnorm(translation_vectors,2,1),'ascend');
            translation_vectors=translation_vectors(:,sort_index);

            % Assign the computed translation_vectors to
            % color_step_translation property
            obj.color_step_translation=translation_vectors;

        end

        function integral_value = IntegrateOverLattice(obj,integrand,index_array)
            % This method computes a 2D integral over a generic lattice
            % by applying the trapezoidal integration rule.
            % INPUTS:
            %   integrand:      a 1xN containing the integrand values
            %   index_array:    a 2xN matrix containing the lattice indexes
            %                   associated to each integrand point
            % OUTPUTS:
            %   integral_value: value of the integral

            % Extract minimum and maximum indexes in the two directions
            n_min=min(index_array(1,:));
            n_max=max(index_array(1,:));
            m_min=min(index_array(2,:));
            m_max=max(index_array(2,:));

            % Create n and m arrays
            n=n_min:n_max;
            m=m_min:m_max;

            % Create n and m grids
            [n,m]=meshgrid(n,m);

            % Extract n and m size
            [Nm,Nn]=size(n);

            % Initialize gridded integrand
            if isa(integrand,'gpuArray')
                gridded_integrand=zeros(Nm,Nn,'gpuArray');
            else
                gridded_integrand=zeros(Nm,Nn);
            end
            
            % Computation of the linear index that maps each index_array
            % couple to its corresponding position on the n, m grids
            linear_index=sub2ind([Nm,Nn],index_array(2,:)-m_min+1,index_array(1,:)-n_min+1);

            % insert the integrand values in the correct grid position
            gridded_integrand(linear_index)=integrand;
            clear integrand linear_index index_array;

            % Compute the integral
            integral_value=trapz(m(:,1),trapz(n(1,:),gridded_integrand,2),1)*abs(det(obj.unit_cell));

        end


    end

end