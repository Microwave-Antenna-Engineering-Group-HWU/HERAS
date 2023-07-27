% This is an example script to simulate a single feed parabolic reflector
% antenna. The results are computed using the HERAS tool and are compared
% against the one obtained with the student edition of GRASP. 

clear all;
close all;
clc;

%% Add paths
addpath(genpath('../../HERAS_tool'));
addpath(genpath('./tempGRASP'));
addpath(genpath('./GRASP_functions'));

%% GRASP parameters
OS='Windows';
GRASPmethod='po';   % 'po' or 'po_plus_ptd'
GRASPSV=1;          % GRASP Student version option
GRASPAR=0;          % GRASP Add Radiation of feed option

% Retrieve path to the unix executable file grasp-analysis on different OSs
% The following lines should automatically identify the path. In case the
% following lines fail to find the path, the user should manually store the
% path to the unix executable file "grasp-analysis" in the variable 
% GRASPpath.
switch OS
    case 'macOS'
        addpath('/Applications');
        script_directory=pwd;
        cd('/Applications')
        list=dir('*GRASP*');
        addpath(genpath(fullfile(list.folder,list.name)));
        GRASPpath=which('grasp-se.idl');
        GRASPpath=strrep(GRASPpath,'grasp-se.idl','grasp-analysis');
        cd(script_directory);
        clear script_directory list
    case 'Windows'
        addpath('C:\Program Files (x86)');
        script_directory=pwd;
        cd('C:\Program Files (x86)')
        list=dir('*TICRA*');
        addpath(genpath(fullfile(list.folder,list.name)));
        GRASPpath=which('grasp-se.idl');
        GRASPpath=strrep(GRASPpath,'grasp-se.idl','grasp-analysis');
        cd(script_directory);
        clear script_directory list        
end

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
T=-12;          % Feed taper

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
feed.setForceFarField(1);

% feed displacement
feed.displaceLocal(displacement_vector);

%% Plot system
f1=figure();
reflector.disp2fig(f1);
feed.disp2fig(f1);
axis tight;

%% Farfield parameters
uv_ThetaPhi=2;
theta_max_deg=90;

% FF parameters
if uv_ThetaPhi == 1
    N_Thetau=51;
    max_Thetau=0.1;
    min_Thetau=-0.1;
    N_Phiv=51;
    min_Phiv=-0.1;
    max_Phiv=0.1;
    Thetau_single=linspace(min_Thetau,max_Thetau,N_Thetau);
    Phiv_single=linspace(min_Phiv,max_Phiv,N_Phiv);
elseif uv_ThetaPhi==2
    N_Thetau=1001;
    max_Thetau=deg2rad(theta_max_deg);
    min_Thetau=deg2rad(-theta_max_deg);
    N_Phiv=3;
    min_Phiv=deg2rad(0);
    max_Phiv=deg2rad(90);
    Thetau_single=linspace(min_Thetau,max_Thetau,N_Thetau);
    Phiv_single=linspace(min_Phiv,max_Phiv,N_Phiv);
end

[Thetau,Phiv]=meshgrid(Thetau_single,Phiv_single);

if uv_ThetaPhi == 1
    cut = FarfieldSphericalGrid();
elseif uv_ThetaPhi == 2
    cut = FarfieldSphericalCut();
end

cut.defineCutsFromSingle(min_Thetau,max_Thetau,N_Thetau,min_Phiv,max_Phiv,N_Phiv);

if uv_ThetaPhi == 1
    UV_grid=cut;
elseif uv_ThetaPhi == 2
    UV_grid=cut.UV_Grid;
end

%% Create coverage
coverage=Coverage(reflector.position(1),reflector.position(2),reflector.position(3),reflector.direction(:,1),reflector.direction(:,2),reflector.direction(:,3));
coverage.DefineCoverageFromUVW(UV_grid.U,UV_grid.V,UV_grid.W);

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

%% Calculate farfield GRASP

% Run GRASP exe file to compute farfields
tic;
[E_CoPol_GRASP,E_XPol_GRASP,time]=CalculateAFR_GRASP(feed,reflector,cut,GRASPmethod,GRASPpath,GRASPSV,GRASPAR);
toc;

% Store GRASP results in Currents2Farfield object
farfield_GRASP = Currents2Farfield();
farfield_GRASP.copolar_component=FF_copol;
eval(sprintf('farfield_GRASP.E_%s=squeeze(E_CoPol_GRASP);',FF_copol));
eval(sprintf('farfield_GRASP.E_%s=squeeze(E_XPol_GRASP);',FF_xpol));
eval(sprintf('farfield_GRASP.D_%s_dB=10*log10(abs(squeeze(E_CoPol_GRASP)).^2/(2*eta_0)*4*pi/feed.Prad);',FF_copol));
eval(sprintf('farfield_GRASP.D_%s_dB=10*log10(abs(squeeze(E_XPol_GRASP)).^2/(2*eta_0)*4*pi/feed.Prad);',FF_xpol));


%% Plot farfields
Theta_single_plot=rad2deg(cut.Theta_single);

% Extract CoPol and XPol components
[E_CoPol_HERAS,E_XPol_HERAS,D_CoPol_dB_HERAS,D_XPol_dB_HERAS]=farfield_HERAS.ExtractCoPolXPolComponents();
[E_CoPol_GRASP,E_XPol_GRASP,D_CoPol_dB_GRASP,D_XPol_dB_GRASP]=farfield_GRASP.ExtractCoPolXPolComponents();

% Phase normalization HERAS
PhaseCompensationFactor=exp(-1i*angle(E_CoPol_HERAS(1,ceil(N_Thetau/2))));
E_CoPol_HERAS = E_CoPol_HERAS *PhaseCompensationFactor;
E_XPol_HERAS  = E_XPol_HERAS  *PhaseCompensationFactor;

% Phase normalization GRASP
PhaseCompensationFactor=exp(-1i*angle(E_CoPol_GRASP(1,ceil(N_Thetau/2))));
E_CoPol_GRASP = E_CoPol_GRASP *PhaseCompensationFactor;
E_XPol_GRASP  = E_XPol_GRASP  *PhaseCompensationFactor;

% HERAS PO plot functions
Phase_CoPol_HERAS = rad2deg(angle(E_CoPol_HERAS));
Phase_XPol_HERAS  = rad2deg(angle(E_XPol_HERAS));

% GRASP PO plot functions
Phase_CoPol_GRASP = rad2deg(angle(E_CoPol_GRASP));
Phase_XPol_GRASP  = rad2deg(angle(E_XPol_GRASP));

%% Figures

% Axis limits
D_max_scale=ceil(max(D_CoPol_dB_GRASP,[],'all')/5)*5;
angle_lim=[min(Theta_single_plot) max(Theta_single_plot)];
Gain_lim=[D_max_scale-120 D_max_scale];
phase_lim=[-180 180];

% for loop over phi angles
for n=1:N_Phiv

    % Gain co-polar component
    fig_amp_copol(n)=figure();
    hold on;
    plot_blue (Theta_single_plot,D_CoPol_dB_GRASP(n,:));
    plot_red  (Theta_single_plot,D_CoPol_dB_HERAS(n,:));
    hold off;
    grid on; grid minor;
    xlim(angle_lim);
    ylim(Gain_lim);
    xlabelBIG('$\theta$ [deg]'); ylabelBIG('Gain [dBi]');
    titleBIG(sprintf('Gain, %s component, $\\Phi = %.1f$ deg',FF_copol,rad2deg(Phiv(n,1))));
    legend('GRASP (PO)','HERAS (PO)');  
    GCAaxisBIG();

    % Gain cross-polar component
    fig_amp_xpol(n)=figure();
    hold on;
    plot_blue (Theta_single_plot,D_XPol_dB_GRASP(n,:));
    plot_red  (Theta_single_plot,D_XPol_dB_HERAS(n,:));
    hold off;
    grid on; grid minor;
    xlim(angle_lim);
    ylim(Gain_lim);
    xlabelBIG('$\theta$ [deg]'); ylabelBIG('Gain [dBi]');
    titleBIG(sprintf('Gain, %s component, $\\Phi = %.1f$ deg',FF_xpol,rad2deg(Phiv(n,1))));
    legend('GRASP (PO)','HERAS (PO)');  
    GCAaxisBIG();

    % Phase co-polar component
    fig_phase_copol(n)=figure();
    hold on;
    plot_blue (Theta_single_plot,Phase_CoPol_GRASP(n,:));
    plot_red  (Theta_single_plot,Phase_CoPol_HERAS(n,:));
    hold off;
    grid on; grid minor;
    xlim(angle_lim);
    ylim(phase_lim);
    xlabelBIG('$\theta$ [deg]'); ylabelBIG('Phase [deg]');
    titleBIG(sprintf('Phase, %s component, $\\Phi = %.1f$ deg',FF_copol,rad2deg(Phiv(n,1))));
    legend('GRASP (PO)','HERAS (PO)');  
    GCAaxisBIG();

    % Phase cross-polar component
    fig_phase_xpol(n)=figure();
    hold on;
    plot_blue (Theta_single_plot,Phase_XPol_GRASP(n,:));
    plot_red  (Theta_single_plot,Phase_XPol_HERAS(n,:));
    hold off;
    grid on; grid minor;
    xlim(angle_lim);
    ylim(phase_lim);
    xlabelBIG('$\theta$ [deg]'); ylabelBIG('Phase [deg]');
    titleBIG(sprintf('Phase, %s component, $\\Phi = %.1f$ deg',FF_xpol,rad2deg(Phiv(n,1))));
    legend('GRASP (PO)','HERAS (PO)');  
    GCAaxisBIG();
end
