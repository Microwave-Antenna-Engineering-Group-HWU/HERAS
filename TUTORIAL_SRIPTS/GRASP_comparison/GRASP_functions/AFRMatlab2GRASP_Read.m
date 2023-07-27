function [E_co,E_cross] = AFRMatlab2GRASP_Read(GRASPTemplate,resultsfilename)
%AFRMATLAB2GRASP_TOR Summary of this function goes here
%   Detailed explanation goes here

LPCP=GRASPTemplate.feed.LPCP;
LHRH=GRASPTemplate.feed.LHRH;
orientation=GRASPTemplate.feed.orientation;
%% 
%% Read results
GRASPTemplate.FarfieldGRASP = Currents2Farfield_old();
if isa(GRASPTemplate.cut,'FarfieldSphericalCut')
    results_filename=fullfile(GRASPTemplate.FolderPath,resultsfilename);
    pattern=readsphcutgrasp(results_filename);
    delete(fullfile(GRASPTemplate.FolderPath, '*.xml'));
    GRASPTemplate.FarfieldGRASP.cut=FarfieldSphericalCut();
    GRASPTemplate.FarfieldGRASP.cut.defineCutsFromSingle(min(pattern{1,2}),max(pattern{1,2}),length(pattern{1,2}),min(pattern{1,3}),max(pattern{1,3}),length(pattern{1,3}))
elseif isa(GRASPTemplate.cut,'FarfieldSphericalGrid')
    results_filename=fullfile(GRASPTemplate.FolderPath,resultsfilename);
    [E_1,E_2,u,v]=readsphgridgrasp(results_filename);
    delete(fullfile(GRASPTemplate.FolderPath, '*.xml'));
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
            phase_co=deg2rad(pattern{1,1}(:,2,i));
            E_co(i,:)=magnitude_co.*exp(1j*phase_co);
            magnitude_cross=sqrt(10.^(pattern{1,1}(:,3,i)/10)*2*eta_0*P_rad/(4*pi));
           phase_cross=deg2rad(pattern{1,1}(:,4,i));
            E_cross(i,:)=magnitude_cross.*exp(1j*phase_cross);
        elseif (isequal(LPCP,'LP') && orientation==pi/2) || (isequal(LPCP,'CP') && isequal(LHRH,'RH'))
            magnitude_co=sqrt(10.^(pattern{1,1}(:,3,i)/10)*2*eta_0*P_rad/(4*pi));
            phase_co=deg2rad(pattern{1,1}(:,4,i));
            E_co(i,:)=magnitude_co.*exp(1j*phase_co);
            magnitude_cross=sqrt(10.^(pattern{1,1}(:,1,i)/10)*2*eta_0*P_rad/(4*pi));
            phase_cross=deg2rad(pattern{1,1}(:,2,i));
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


end

