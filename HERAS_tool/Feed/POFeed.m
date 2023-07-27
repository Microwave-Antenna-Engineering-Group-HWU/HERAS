classdef (Abstract) POFeed < POobject
    % POFeed describes a general feed. It defines the interfaces to be
    %   implemented by any feed to integrate seemlessly with the other parts
    %   of the PO tool, and provide basic methods common to all feeds
    %   It inherits from the abstract class POobject
    %   Author : Louis Dufour
    %   Revisions : 0.1 - 25/10/2018
    %
    % POFeed methods (Abstract)
    %    illuminate : given a set of points, give the H fields obtained by
    %          illumination from the feed. All points and vectors in feed
    %          CS
    %
    % POFeed methods
    %   displayFeed : display the feed either in the figure given in
    %   parameters, or by creating a new figure, with the CS. The CS scale
    %   is given by POobject.CSscale
    %
    % See POobject

    properties (Constant)
        theta1 = 0.05;
        mu_0=4*pi*1e-7;      % Vacuum permeability
        epsilon_0=8.854e-12; % Vacuum permittivity
    end
    
    properties
        k_0, eta_0, phase_yn, freq;
        Prad;
    end
    
    methods (Abstract)
        [E_x_s,E_y_s,E_z_s,H_x_s,H_y_s,H_z_s] = illuminate(obj,x_s,y_s,z_s);
        disp2fig(obj,ax, CSref);
    end
        
end