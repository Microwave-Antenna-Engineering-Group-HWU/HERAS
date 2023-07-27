function AFRMatlab2GRASP_TCI(GRASPTemplate,Array,GRASPSV,GRASPAR)
%AFRMATLAB2GRASP_TOR Summary of this function goes here
%   Detailed explanation goes here

f=GRASPTemplate.reflector.surface.f;
x_c=GRASPTemplate.reflector.surface.x_c;
z_f=-f*(-1+x_c^2/(4*f^2));

%% Set TCI tasks
data_tci{(length(Array)-1)*4+4,:}=newline;
for i=1:length(Array)
    data_tci{(i-1)*4+1,:}=sprintf('COMMAND OBJECT single_feed_coor set ( origin : struct(x: %.12g m, y: %.12g m, z:  %.12g m), class : coor_sys) cmd_%d ',Array(i).position(1)+x_c,Array(i).position(2),f-z_f+Array(i).position(3),(i-1)*4+1);
    if isa(GRASPTemplate.cut,'FarfieldSphericalCut')
        data_tci{(i-1)*4+2,:}=sprintf('COMMAND OBJECT single_cut set ( file_name : single_cut%d.cut, class : spherical_cut)  cmd_%d',i,(i-1)*4+2);
        if GRASPSV==1
            data_tci{(i-1)*4+3,:}=sprintf('COMMAND OBJECT single_po get_currents ( source : sequence(ref(single_feed)))  cmd_%d',(i-1)*4+3);
        else
            data_tci{(i-1)*4+3,:}=sprintf('COMMAND OBJECT single_po get_currents ( source : sequence(ref(single_feed)), auto_convergence_of_po : on, convergence_on_output_grid : ref(spherical_cut))  cmd_%d',(i-1)*4+3);
        end
        if GRASPAR == 1
            data_tci{(i-1)*4+4,:}=sprintf('COMMAND OBJECT single_cut get_field ( source : sequence(ref(single_feed),  ref(single_po))) cmd_%d',(i-1)*4+4);
        else
            data_tci{(i-1)*4+4,:}=sprintf('COMMAND OBJECT single_cut get_field ( source : sequence(ref(single_po))) cmd_%d',(i-1)*4+4);
        end
    elseif isa(GRASPTemplate.cut,'FarfieldSphericalGrid')
        data_tci{(i-1)*4+2,:}=sprintf('COMMAND OBJECT spherical_grid set ( file_name : spherical_grid%d.grd, class : spherical_grid)  cmd_%d',i,(i-1)*4+2);
        if GRASPSV==1
            data_tci{(i-1)*4+3,:}=sprintf('COMMAND OBJECT single_po get_currents ( source : sequence(ref(single_feed)))  cmd_%d',(i-1)*4+3);
        else
            data_tci{(i-1)*4+3,:}=sprintf('COMMAND OBJECT single_po get_currents ( source : sequence(ref(single_feed)), auto_convergence_of_po : on, convergence_on_output_grid : ref(spherical_grid))  cmd_%d',(i-1)*4+3);
        end
        if GRASPAR == 1
            data_tci{(i-1)*4+4,:}=sprintf('COMMAND OBJECT spherical_grid get_field ( source : sequence(ref(single_feed),  ref(single_po))) cmd_%d',(i-1)*4+4);
        else
            data_tci{(i-1)*4+4,:}=sprintf('COMMAND OBJECT spherical_grid get_field ( source : sequence(ref(single_po))) cmd_%d',(i-1)*4+4);
        end
    end
end
data_tci{(i-1)*4+5,:}=sprintf('QUIT ');
data_tci{(i-1)*4+6,:}=newline;

%% Write new TCI file
filename=fullfile(GRASPTemplate.FolderPath, GRASPTemplate.TCIfilename);
delete(filename);
newfilename=fullfile(GRASPTemplate.FolderPath, GRASPTemplate.TCIfilename);
fileID = fopen(newfilename,'w');
for i=1:size(data_tci,1)
    fprintf(fileID,[data_tci{i,:} '\r\n']);
end
fclose(fileID);

end

