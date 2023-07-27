classdef Lattice < matlab.mixin.Copyable
    % Lattice is a class implementing a general N dimensional lattice. This
    % class supports multiple colours defined as sub-lattices of the
    % original lattice.
    %
    % Author : Francesco Lisi
    % Revisions : 0.1.0 - 28/02/2023

    properties
        N_dimensions;   % N_dimensions defines the lattice dimensions 
        unit_cell;      % lattice unit cell (N_dimensions x N_dimensions matrix)
        N_colors;
        M_colors_matrix;
        color_step_translation;
    end
        
    methods (Abstract)
        SetMultipleColors(obj,N_colors);
    end

    methods

        function obj = Lattice(unit_cell)
            % Constructor

            % Store unit cell size
            unit_cell_size=size(unit_cell);

            % check if unit_cell is a square matrix
            if (length(unit_cell_size)==2)&&(unit_cell_size(1)==unit_cell_size(2))
                % Assign properties to the object
                obj.N_dimensions=unit_cell_size(1);
                obj.unit_cell=unit_cell;
            else
                error('wrong input size: input should be a square matrix.')
            end
        end

        function points=getLatticePoints(obj,n)
            % This method computes the lattice points in cartesian
            % coordinates from the index vector n.

            % Check that n is a vector of integers
            if sum(abs(mod(n,1)),'all')==0
                % Compute points
                points=obj.unit_cell*n;
            else
                error('Wrong input: the input should be an array of integers');
            end
        end
    
        function [points,varargout] = GeneratePointsInsideContour(obj,contour)
            % This function generates all the points that belong to the
            % lattice and are inside the region defined by contour. points
            % is the location of the point in the cartesian grid, while
            % step_array corresponds to points=unit_cell*step_array

            if contour.N_dimensions~=obj.N_dimensions
                error('Incompatible input contour: the contour dimension must be the same as the lattice one.')
            else

                % Define color unit cell
                unit_cell_color=obj.unit_cell*obj.M_colors_matrix;                

                % minimum distance between adjacent points
                minimum_interelement_distance=min(vecnorm(unit_cell_color,2,1));
    
                % maximum number of allowed steps. This parameter ensures that
                % we generate all the points inside the contour, but we do not
                % create a step_array with too many elements outside it, since
                % they will be removed.
                maximum_steps=2*ceil(contour.furthest_point/minimum_interelement_distance)+max(sum(abs(obj.M_colors_matrix),1));
    
                %% Compute lattice points
                % step_array is a matrix whose columns correspond to the number 
                % of steps to take in each principal direction defined by 
                % unit_cell.
    
                % First, we generate all the possible combination of steps,
                % where the maximum number of steps along each direction
                % corresponds to maximum_steps
                step_array_color=zeros(obj.N_dimensions,(2*maximum_steps+1)^obj.N_dimensions);
                for n=1:obj.N_dimensions
                    step_array_color(n,:)=repmat(kron((-maximum_steps:maximum_steps),ones(1,(2*maximum_steps+1)^(obj.N_dimensions-n))),[1,(2*maximum_steps+1)^(n-1)]);
                end
    
                % Then, we remove all the columns with a total number of steps
                % greater than maximum_steps
                step_array_color=step_array_color(:,sum(abs(step_array_color))<=maximum_steps);
    
                % Computation of the corresponding points in the cartesian grid
                points_color=unit_cell_color*step_array_color;

                % If N_colors==1, points and step_array are arrays. If
                % N_colors>1 then points and step_array are cell arrays,
                % where each cell element corresponds to the points and
                % step_array associated to that color.
                if obj.N_colors==1
    
                    %% Remove points outside the contour
                    % Obtain points inside the contour
                    inside_contour_mask=contour.IsInside(points_color);
        
                    % remove elements outside contour
                    points=points_color(:,inside_contour_mask);
                    step_array=step_array_color(:,inside_contour_mask);
                    
                else

                    % Initialize points and step_array cell arrays. each
                    % cell array element corresponds to one color
                    points=cell(obj.N_colors,1);
                    step_array=cell(obj.N_colors,1);

                    for c=1:obj.N_colors

                        % Compute step_array and points for each color
                        step_array{c}=obj.M_colors_matrix*step_array_color+obj.color_step_translation(:,c);
                        points{c}=obj.unit_cell*step_array{c};

                        %% Remove points outside the contour
                        % Obtain points inside the contour
                        inside_contour_mask=contour.IsInside(points{c});
            
                        % remove elements outside contour
                        step_array{c}=step_array{c}(:,inside_contour_mask);
                        points{c}=points{c}(:,inside_contour_mask);

                    end

                end

                % the second output correponds to step_array where the elements
                % outside the contour are removed
                if nargout==2
                    varargout{1}=step_array;
                end

            end
        end

    end

end
