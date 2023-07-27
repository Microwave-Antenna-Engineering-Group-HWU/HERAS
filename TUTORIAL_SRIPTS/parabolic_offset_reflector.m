% This is an example script to simulate a single feed parabolic reflector
% antenna. 

clear all;
close all;
clc;

%% Add paths
addpath(genpath('../HERAS_tool'));

%% Electrical parameters
f_0=17.7e9;                 % Frequency
epsilon_0=8.854e-12;        % Vacuum permittivity
mu_0=4*pi*1e-7;             % Vacuum permeability
c_0=physconst('LightSpeed');% speed of light
Y_0=sqrt(epsilon_0/mu_0);   % Admittance of free space
eta_0=1/Y_0;                % Impedance of free space
lambda_0=c_0/f_0;           % wavelength
k_0=2*pi/lambda_0;          % wavenumber

%% Main reflector
D=2.75;     % Reflector rim diameter
F=2*D;      % Focal length
x_c=3.1;    % Offset height
h=x_c-D/2;  % Offset hight (clearance)

% Create reflector object
reflector = CircleReflector(D,0,0,0);

% Create reflector surface
surface = ParabSurface(F, x_c, h);

% Assign surface to reflector
reflector.setSurface(surface);

%% Mesh

% Compute number of rho and phi elements for the mesh
theta_max=deg2rad(90);                      % Maximum angle for convergence
aux=1.09*pi*(D/lambda_0)*sin(theta_max)+10;
Nrho=round(aux/2.4);
Nphi=round(aux);

% Add more points for convergence
Nrho=ceil(Nrho*1.25); 
Nphi=ceil(Nphi*1.25);

% Create mesh
mesh=PolarMeshGRASP(Nrho, Nphi);

% Assign mesh to reflector
reflector.setMesh(mesh);

% Mesh and apply surface to reflector
reflector.meshReflector();
reflector.depthReflector();

%% Feed
phase_yn=1;     % always set to 1
LPCP='CP';      % 'LP' for linear polarization; 'CP' for circular
LHRH='LH';      % Sense of rotation for 'CP' polarization
orientation=0;  % Feed orientation expressed as the anticlockwise rotation angle with respect to x direction
T=-12;          % Feeds taper

% Feed offset respect to the focus
displacement_vector=[0;0;0];

% The follwing lines of code determine the co-polar and cross-polar
% components of the fields based on the previous parameters
switch LPCP
    case 'LP'
        switch orientation
            case 0
                FF_copol='V';
                FF_xpol ='H';

            case pi/2
                FF_copol='H';
                FF_xpol ='V';
        end

    case 'CP'
        LPCP_2='CP';
        switch LHRH
            case 'LH'
                FF_copol='RHCP';
                FF_xpol ='LHCP';
            case 'RH'
                FF_copol='LHCP';
                FF_xpol ='RHCP';
        end
end

% Distance along z direction between the focal point and the center of the
% reflector
z_f=-F*(-1+x_c^2/(4*F^2));

% Feed CS direction and position with respect to Global CS
feed_ez = -[-x_c,0, z_f]';feed_ey = [0,-1,0]'; feed_ex = cross(feed_ey, feed_ez);
feed_pos = [-x_c,0,z_f]';

% Ideal gaussian source initialization
feed=IdealGaussianSourceFeed(feed_pos(1), feed_pos(2), feed_pos(3),feed_ex, feed_ey, feed_ez);

% Compute angle at the lower rim of the reflector
theta_feed=reflector.getTaperAngle(feed,'low');

% Set feed parameteres
feed.setParameters(T,theta_feed,LPCP,LHRH,orientation,phase_yn,k_0,eta_0);

% feed.forceFarField == 1 -> Feed FF expression (Default);
% feed.forceFarField == 0 -> Feed NF expression;
feed.setForceFarField(0);

% feed displacement
feed.displaceLocal(displacement_vector);

%% Plot system
f1=figure();
reflector.disp2fig(f1);
feed.disp2fig(f1);
axis tight;

%% Create coverage

% Number of theta and phi points
N_theta=1001;
N_phi=3;

% extremes of the intervals
max_theta = deg2rad(90);
min_theta = -max_theta;
max_phi   = deg2rad(90);
min_phi   = 0;

% Theta and Phi mesh creation
Theta   = linspace(min_theta,max_theta,N_theta);
Phi     = linspace(min_phi,max_phi,N_phi);
[Theta,Phi]=meshgrid(Theta,Phi);

% Create coverage
coverage=Coverage(reflector.position(1),reflector.position(2),reflector.position(3),reflector.direction(:,1),reflector.direction(:,2),reflector.direction(:,3));
coverage.DefineCoverageFromThetaPhi(Theta,Phi);

%% Calculate farfield

% Create surface currents object
currents = POSurfaceCurrents(reflector,feed);

% Compute surface currents
currents.calculateCurrents();

% Create farfield object
farfield_HERAS = Currents2Farfield();

% Compute farfield
tic;
farfield_HERAS.calculateFarField(currents,coverage,FF_copol);
toc;

%% Plot farfields
Theta_single_plot=rad2deg(Theta(1,:));

% Extract CoPol and XPol components
[E_CoPol_HERAS,E_XPol_HERAS,D_CoPol_dB_HERAS,D_XPol_dB_HERAS]=farfield_HERAS.ExtractCoPolXPolComponents();

% Phase normalization HERAS
PhaseCompensationFactor=exp(-1i*angle(E_CoPol_HERAS(1,ceil(N_theta/2))));
E_CoPol_HERAS = E_CoPol_HERAS *PhaseCompensationFactor;
E_XPol_HERAS  = E_XPol_HERAS  *PhaseCompensationFactor;

% HERAS PO plot functions
Phase_CoPol_HERAS = rad2deg(angle(E_CoPol_HERAS));
Phase_XPol_HERAS  = rad2deg(angle(E_XPol_HERAS));

%% Figures

% Axis limits
D_max_scale=ceil(max(D_CoPol_dB_HERAS,[],'all')/5)*5;
angle_lim=[min(Theta_single_plot) max(Theta_single_plot)];
Gain_lim=[D_max_scale-120 D_max_scale];
phase_lim=[-180 180];

% for loop over phi angles
for n=1:N_phi

    % Gain co-polar component
    fig_amp(n)=figure();
    hold on;
    plot_red  (Theta_single_plot,D_CoPol_dB_HERAS(n,:));
    plot_blue (Theta_single_plot,D_XPol_dB_HERAS(n,:));
    hold off;
    grid on; grid minor;
    xlim(angle_lim);
    ylim(Gain_lim);
    xlabelBIG('$\theta$ [deg]'); ylabelBIG('Gain [dBi]');
    titleBIG(sprintf('Gain, $\\phi = %.1f$ deg',rad2deg(Phi(n,1))));
    legend([string(FF_copol), string(FF_xpol)])
    GCAaxisBIG();

    % Phase co-polar component
    fig_phase_copol(n)=figure();
    hold on;
    plot_red  (Theta_single_plot,Phase_CoPol_HERAS(n,:));
    hold off;
    grid on; grid minor;
    xlim(angle_lim);
    ylim(phase_lim);
    xlabelBIG('$\theta$ [deg]'); ylabelBIG('Phase [deg]');
    titleBIG(sprintf('Phase, %s component, $\\phi = %.1f$ deg',FF_copol,rad2deg(Phi(n,1))));
    GCAaxisBIG();

    % Phase cross-polar component
    fig_phase_xpol(n)=figure();
    hold on;
    plot_blue  (Theta_single_plot,Phase_XPol_HERAS(n,:));
    hold off;
    grid on; grid minor;
    xlim(angle_lim);
    ylim(phase_lim);
    xlabelBIG('$\theta$ [deg]'); ylabelBIG('Phase [deg]');
    titleBIG(sprintf('Phase, %s component, $\\phi = %.1f$ deg',FF_xpol,rad2deg(Phi(n,1))));
    GCAaxisBIG();
end
