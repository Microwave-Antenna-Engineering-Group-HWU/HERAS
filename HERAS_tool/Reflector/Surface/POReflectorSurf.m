classdef (Abstract) POReflectorSurf < matlab.mixin.Copyable
    % POReflectorSurf is the abstract class (interface) for all
    % implementations (subclases) of reflector surfaces to follow
    %
    %   Author : Louis Dufour
    %   Revisions : 0.1.0 - 25/10/2018
    %
    % methods (Abstract)
    %   getSurf     - for a set of x,y returns the z. Vectorized.
    %   getNorm     - for a set of x,y returns the normals and Jacobian(?). Vectorized.
    %   getLines    - returns the lines that are forced as edges of the
    %                mesh. Note that as of v0.1.0 the meshing toolbox used
    %                does not deal well with it and creates really small
    %                mesh cells at the edge, which considerably increase
    %                the computation time without adding any significant
    %                advantage at this stage (this stage = limited
    %                performances optimisation)
    
    properties
    end
    
    methods (Abstract)
        z = getSurf(obj, x, y)
        [n_x, n_y, n_z, N] = getNorm(obj, x, y)
        lines = getLines(obj,x,y);
    end
    
end
