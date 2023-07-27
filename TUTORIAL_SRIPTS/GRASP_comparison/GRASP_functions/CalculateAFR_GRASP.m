function [Array_farfields_CO,Array_farfields_XP,time]=CalculateAFR_GRASP(Array,reflector,cut,GRASPmethod,GRASPpath,GRASPSV,GRASPAR)
%CALCULATEAFR_GRASP Summary of this function goes here
%   Detailed explanation goes here

if  isa(cut,'FarfieldSphericalGrid')
    %% Correct U-V grid to GRASP definition
    u_single=unique(cut.U(~isnan(cut.U(:))),'sorted');
    v_single=unique(cut.V(~isnan(cut.V(:))),'sorted');

    [u_new,v_new]=meshgrid(u_single,v_single); %Rectangular grid for GRASP

    %Position corrections between GRASP rectangular Grid and Matlab hexagonal
    %grid
    pos=nan(size(cut.U));
    for i=1:length(cut.U(:))
        if ~isnan(cut.U(i))
            pos(i)=find(abs(cut.U(i)-u_new)<1e-5 & abs(cut.V(i)-v_new)<1e-5);  
        end
    end

    cutGRASP=FarfieldSphericalGrid();
    cutGRASP.defineCutsFromSingle(min(u_single),max(u_single),length(u_single),min(v_single),max(v_single),length(v_single));

elseif isa(cut,'FarfieldSphericalCut')
    cutGRASP=cut.copy();
end

%% Write GRASP files and run simulation
t1=datetime;
GRASPTemplate=SingleReflectorGRASPTemplate(reflector,Array(1),cutGRASP,GRASPmethod);
GRASPTemplate.FolderPath=fullfile(pwd,'tempGRASP');
AFRMatlab2GRASP_TOR(GRASPTemplate);
AFRMatlab2GRASP_TCI(GRASPTemplate,Array,GRASPSV,GRASPAR);
t2=datetime;
time.WriteGRASP=t2-t1;
t1=datetime;
AFRMatlab2GRASP_Run(GRASPTemplate,GRASPpath);
t2=datetime;
time.RunGRASP=t2-t1;

if  isa(cut,'FarfieldSphericalGrid')
    
    E_CO_GRASP=nan(length(Array),size(cutGRASP.U,1),size(cutGRASP.U,2));
    E_XP_GRASP=nan(length(Array),size(cutGRASP.U,1),size(cutGRASP.U,2));

    t1=datetime;
    for i=1:length(Array)
        resultsfilename=['spherical_grid' num2str(i) '.grd'];
        [E_CO_GRASP(i,:,:),E_XP_GRASP(i,:,:)]=AFRMatlab2GRASP_Read(GRASPTemplate,resultsfilename);
        delete(fullfile(GRASPTemplate.FolderPath, resultsfilename));
    end
    t2=datetime;
    time.ReadGRASP=t2-t1;

elseif isa(cut,'FarfieldSphericalCut')
    
    E_CO_GRASP=nan(length(Array),size(cutGRASP.Theta,1),size(cutGRASP.Theta,2));
    E_XP_GRASP=nan(length(Array),size(cutGRASP.Theta,1),size(cutGRASP.Theta,2));

    t1=datetime;
    for i=1:length(Array)
        resultsfilename=['single_cut' num2str(i) '.cut'];
        [E_CO_GRASP(i,:,:),E_XP_GRASP(i,:,:)]=AFRMatlab2GRASP_Read(GRASPTemplate,resultsfilename);
        delete(fullfile(GRASPTemplate.FolderPath, resultsfilename));
    end
    t2=datetime;
    time.ReadGRASP=t2-t1;
end

%% Convert to final grid

if  isa(cut,'FarfieldSphericalGrid')

    Array_farfields_CO=nan(length(Array),length(find(~isnan(cut.U(:)))));
    Array_farfields_XP=nan(length(Array),length(find(~isnan(cut.U(:)))));

    for i=1:length(Array)
        aux=nan(size(cut.U));
        aux(~isnan(cut.U))=E_CO_GRASP(i,pos(~isnan(cut.U(:))));
        Array_farfields_CO(i,:)=aux(~isnan(cut.U(:)));
        aux=nan(size(cut.U));
        aux(~isnan(cut.U))=E_XP_GRASP(i,pos(~isnan(cut.U(:))));
        Array_farfields_XP(i,:)=aux(~isnan(cut.U(:)));
    end

elseif isa(cut,'FarfieldSphericalCut')

    Array_farfields_CO=E_CO_GRASP;
    Array_farfields_XP=E_XP_GRASP;
    
end


end

