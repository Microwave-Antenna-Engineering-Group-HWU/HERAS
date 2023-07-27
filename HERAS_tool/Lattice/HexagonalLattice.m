classdef HexagonalLattice < Lattice2d
    % This class defines an hexagonal lattice with distance between two
    % adjacent cells equal to intercell_distance

    properties
        intercell_distance;
    end
 
    methods

        function obj = HexagonalLattice(intercell_distance)
            % Constructor
            unit_cell = intercell_distance*[cos(-pi/6), cos(+pi/6); ...
                                            sin(-pi/6), sin(+pi/6)];
            obj@Lattice2d(unit_cell);
            obj.intercell_distance=intercell_distance;
        end

    end

end