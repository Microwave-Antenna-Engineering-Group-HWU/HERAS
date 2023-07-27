function GRASPTemplate = SingleReflectorMatlab2GRASP(GRASPTemplate)
%SINGLEREFLECTORMATLAB2GRASP Summary of this function goes here
%   Detailed explanation goes here

%% Obtain variables from classes
freq=GRASPTemplate.feed.freq;
f=GRASPTemplate.reflector.surface.f;
x_c=GRASPTemplate.reflector.surface.x_c;
a1=GRASPTemplate.reflector.rim.a1;
a2=GRASPTemplate.reflector.rim.a2;
feed_position=GRASPTemplate.feed.position;
z_f=-f*(-1+x_c^2/(4*f^2));
feed_direction=GRASPTemplate.feed.direction;

LPCP=GRASPTemplate.feed.LPCP;
LHRH=GRASPTemplate.feed.LHRH;
orientation=GRASPTemplate.feed.orientation;

%% Read TOR file
filename=[GRASPTemplate.FolderPath GRASPTemplate.TORfilename];

fileID = fopen(filename);
line=fgetl(fileID);
i=1;
while(line~=-1)
    data{i,:}=line;
    i=i+1;
    line=fgetl(fileID);
end
fclose(fileID);

%% Modify file parameters
index=find(contains(data,'single_frequencies  frequency  '));
data{index+2,:}=sprintf('  frequency_list   : sequence(%.12g GHz)',freq/1e9);

index=find(contains(data,'single_surface  paraboloid  '));
data{index+2,:}=sprintf('  focal_length     : %.12g m',f);

index=find(contains(data,'single_rim  elliptical_rim  '));
data{index+2,:}=sprintf('  centre           : struct(x: %.12g m, y: %.12g m),',x_c+GRASPTemplate.reflector.rim.x0,GRASPTemplate.reflector.rim.y0);
data{index+3,:}=sprintf('  half_axis        : struct(x: %.12g m, y: %.12g m)',a1,a2);

index=find(contains(data,'single_feed_coor  coor_sys  '));
data{index+2,:}=sprintf('  origin           : struct(x: %.12g m, y: %.12g m, z: %.12g m),',feed_position(1)+x_c,feed_position(2),f-z_f+feed_position(3));
data{index+3,:}=sprintf('  x_axis           : struct(x: %.12g, y: %.12g, z: %.12g),',feed_direction(1,1),feed_direction(1,2),feed_direction(1,3));
data{index+4,:}=sprintf('  y_axis           : struct(x: %.12g, y: %.12g, z: %.12g),',feed_direction(2,1),feed_direction(2,2),feed_direction(2,3));

index=find(contains(data,'single_feed  gaussian_beam_pattern  '));
data{index+4,:}=sprintf('  taper_angle      : %.12g,',rad2deg(GRASPTemplate.feed.tapering_angle));
data{index+5,:}=sprintf('  taper            : %.12g,',GRASPTemplate.feed.tapering_dB);
if GRASPTemplate.feed.forceFarField
    data{index+6,:}=sprintf('  far_forced       : on,');
else
    data{index+6,:}=sprintf('  far_forced       : off,');
end
switch LPCP
    case 'LP'
        switch orientation
            case 0
                data{index+7,:}=sprintf('  polarisation     : linear_x');
            case pi/2
                data{index+7,:}=sprintf('  polarisation     : linear_y');
            otherwise
                warning('Polarization angle not enable with GRASP');
        end
    case 'CP'
        switch LHRH
            case 'RH'
                data{index+7,:}=sprintf('  polarisation     : rhc');
            case 'LH'
                data{index+7,:}=sprintf('  polarisation     : lhc');
        end
end

index=find(contains(data,'single_cut_coor  coor_sys  '));
data{index+2,:}=sprintf('  origin           : struct(x: %.12g m, y: %.12g m, z: %.12g m),',x_c,0,f-z_f);

index=find(contains(data,'single_po  po_single_face_scatterer  '));
data{index+4,:}=sprintf('  method           : %s',GRASPTemplate.GRASPmethod);

if isa(GRASPTemplate.cut,'FarfieldSphericalGrid')
    u_single=GRASPTemplate.cut.u_single;
    v_single=GRASPTemplate.cut.v_single;
    index=find(contains(data,'spherical_grid  spherical_grid  '));
    data{index+3,:}=sprintf('  x_range          : struct(start: %E, end: %E, np: %d),',min(u_single),max(u_single),length(u_single));
    data{index+4,:}=sprintf('  y_range          : struct(start: %E, end: %E, np: %d),',min(v_single),max(v_single),length(v_single));
    
    switch LPCP
        case 'LP'
            data{index+5,:}=sprintf('  polarisation     : linear,');
        case 'CP'
            data{index+5,:}=sprintf('  polarisation     : circular,');
    end
elseif isa(GRASPTemplate.cut,'FarfieldSphericalCut')
    Theta_single=GRASPTemplate.cut.Theta_single;
    Phi_single=GRASPTemplate.cut.Phi_single;
    index=find(contains(data,'single_cut  spherical_cut  '));
    data{index+3,:}=sprintf('  theta_range      : struct(start: %.12g, end: %.12g, np: %d),',rad2deg(min(Theta_single)),rad2deg(max(Theta_single)),length(Theta_single));
    data{index+4,:}=sprintf('  phi_range        : struct(start: %.12g, end: %.12g, np: %d),',rad2deg(min(Phi_single)),rad2deg(max(Phi_single)),length(Phi_single));
    switch LPCP
        case 'LP'
            data{index+5,:}=sprintf('  polarisation     : linear,');
        case 'CP'
            data{index+5,:}=sprintf('  polarisation     : circular,');
    end
end

%% Write new TOR file
filename=[GRASPTemplate.FolderPath GRASPTemplate.TORfilename];
delete(filename);
newfilename=[GRASPTemplate.FolderPath GRASPTemplate.TORfilename];
fileID = fopen(newfilename,'w');
for i=1:size(data,1)
    fprintf(fileID,[data{i,:} '\r\n']);
end
fclose(fileID);

%% Set TCI tasks
data_tci{1,:}=sprintf('COMMAND OBJECT single_po get_currents ( source : sequence(ref(single_feed)))  &');
data_tci{2,:}=sprintf('single_po_cmd ');
if isa(GRASPTemplate.cut,'FarfieldSphericalCut')
    data_tci{3,:}=sprintf('COMMAND OBJECT single_cut get_field ( source : sequence(ref(single_po)))  &');
elseif isa(GRASPTemplate.cut,'FarfieldSphericalGrid')
    data_tci{3,:}=sprintf('COMMAND OBJECT spherical_grid get_field ( source : sequence(ref(single_po)))  &');
end
data_tci{4,:}=sprintf('single_get_field_cmd ');
data_tci{5,:}=sprintf('QUIT ');
data_tci{6,:}=newline;

%% Write new TCI file
filename=[GRASPTemplate.FolderPath GRASPTemplate.TCIfilename];
delete(filename);
newfilename=[GRASPTemplate.FolderPath GRASPTemplate.TCIfilename];
fileID = fopen(newfilename,'w');
for i=1:size(data_tci,1)
    fprintf(fileID,[data_tci{i,:} '\r\n']);
end
fclose(fileID);

%% Launch batch command
projectdirectory=GRASPTemplate.FolderPath;
work_directory=pwd;
cd(projectdirectory);
[status,cmdout]=system('"C:/Program Files (x86)/TICRA/GRASP-SE-10.3.0/bin/grasp-analysis" batch.gxp main.out main.log');
% [status,cmdout]=system('grasp-analysis batch.gxp main.out main.log');
cd(work_directory);

%% Read results
GRASPTemplate.FarfieldGRASP = Currents2Farfield();
if isa(GRASPTemplate.cut,'FarfieldSphericalCut')
    results_filename=[GRASPTemplate.FolderPath 'single_cut.cut'];
    pattern=readsphcutgrasp(results_filename);
    delete([GRASPTemplate.FolderPath '*.xml']);
    GRASPTemplate.FarfieldGRASP.cut=FarfieldSphericalCut();
    GRASPTemplate.FarfieldGRASP.cut.defineCutsFromSingle(min(pattern{1,2}),max(pattern{1,2}),length(pattern{1,2}),min(pattern{1,3}),max(pattern{1,3}),length(pattern{1,3}))
elseif isa(GRASPTemplate.cut,'FarfieldSphericalGrid')
    results_filename=[GRASPTemplate.FolderPath 'spherical_grid.grd'];
    [E_1,E_2,u,v]=readsphgridgrasp(results_filename);
    delete([GRASPTemplate.FolderPath '*.xml']);
    GRASPTemplate.FarfieldGRASP.cut=FarfieldSphericalGrid();
    GRASPTemplate.FarfieldGRASP.cut.defineCutsFromSingle(min(u(1,:)),max(u(1,:)),length(u(1,:)),min(v(:,1)'),max(v(:,1)'),length(v(:,1)'))
end

%% Make correction
P_rad=GRASPTemplate.feed.Prad;
eta_0=GRASPTemplate.feed.eta_0;
if isa(GRASPTemplate.cut,'FarfieldSphericalCut')
    E_co=zeros(length(pattern{1,3}),length(pattern{1,2}));
    E_cross=zeros(length(pattern{1,3}),length(pattern{1,2}));
    for i=1:size(pattern{1,1},3)
        if (isequal(LPCP,'LP') && orientation==0) || (isequal(LPCP,'CP') && isequal(LHRH,'LH'))
            magnitude_co=sqrt(10.^(pattern{1,1}(:,1,i)/10)*2*eta_0*P_rad/(4*pi));
            phase_co=wrapToPi(deg2rad(pattern{1,1}(:,2,i))+pi);
            E_co(i,:)=magnitude_co.*exp(1j*phase_co);
            magnitude_cross=sqrt(10.^(pattern{1,1}(:,3,i)/10)*2*eta_0*P_rad/(4*pi));
            phase_cross=wrapToPi(deg2rad(pattern{1,1}(:,4,i))+pi);
            E_cross(i,:)=magnitude_cross.*exp(1j*phase_cross);
        elseif (isequal(LPCP,'LP') && orientation==pi/2) || (isequal(LPCP,'CP') && isequal(LHRH,'RH'))
            magnitude_co=sqrt(10.^(pattern{1,1}(:,3,i)/10)*2*eta_0*P_rad/(4*pi));
            phase_co=wrapToPi(deg2rad(pattern{1,1}(:,4,i))+pi);
            E_co(i,:)=magnitude_co.*exp(1j*phase_co);
            magnitude_cross=sqrt(10.^(pattern{1,1}(:,1,i)/10)*2*eta_0*P_rad/(4*pi));
            phase_cross=wrapToPi(deg2rad(pattern{1,1}(:,2,i))+pi);
            E_cross(i,:)=magnitude_cross.*exp(1j*phase_cross);
        end

    end
elseif isa(GRASPTemplate.cut,'FarfieldSphericalGrid')   
    E_1=E_1*sqrt(2*eta_0*P_rad/(4*pi));
    E_2=E_2*sqrt(2*eta_0*P_rad/(4*pi));
    if (isequal(LPCP,'LP') && orientation==0) || (isequal(LPCP,'CP') && isequal(LHRH,'LH'))
        GRASPTemplate.FarfieldGRASP.E_RHCP=E_1;
        E_co=E_1;
        E_cross=E_2;
    elseif (isequal(LPCP,'LP') && orientation==pi/2) || (isequal(LPCP,'CP') && isequal(LHRH,'RH'))
        E_co=E_2;
        E_cross=E_1;
    end
end

%% Save field

switch LPCP
    case 'LP'
        switch orientation
            case 0
                GRASPTemplate.FarfieldGRASP.E_H=E_co;
                GRASPTemplate.FarfieldGRASP.E_V=E_cross;
            case pi/2
                GRASPTemplate.FarfieldGRASP.E_V=E_co;
                GRASPTemplate.FarfieldGRASP.E_H=E_cross;
            otherwise
                warning('Polarization angle not enable with GRASP');
        end
    case 'CP'
        switch LHRH
            case 'RH'
                GRASPTemplate.FarfieldGRASP.E_LHCP=E_co;
                GRASPTemplate.FarfieldGRASP.E_RHCP=E_cross;
            case 'LH'
                GRASPTemplate.FarfieldGRASP.E_RHCP=E_co;
                GRASPTemplate.FarfieldGRASP.E_LHCP=E_cross;
        end
end

end

