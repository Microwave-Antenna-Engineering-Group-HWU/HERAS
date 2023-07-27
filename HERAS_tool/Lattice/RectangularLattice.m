classdef RectangularLattice < Lattice2d
    % This class defines a rectangular lattice with distance between two
    % adjacent cells equal to intercell_distance

    properties
        intercell_distance_x;
        intercell_distance_y;
    end

    methods
        
        function obj = RectangularLattice(varargin)
            % Constructor
            switch nargin
                case 1
                    intercell_distance_x=varargin{1};
                    intercell_distance_y=varargin{1};
                case 2
                    intercell_distance_x=varargin{1};
                    intercell_distance_y=varargin{2};
                otherwise
                    error('Too many input arguments.');
            end                    
            unit_cell = [intercell_distance_x, 0; ...
                         0, intercell_distance_y];
            obj@Lattice2d(unit_cell);
            obj.intercell_distance_x = intercell_distance_x;
            obj.intercell_distance_y = intercell_distance_y;

        end

    end

end