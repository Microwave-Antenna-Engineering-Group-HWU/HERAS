function AFRMatlab2GRASP_Run(GRASPTemplate,GRASPpath)
%AFRMATLAB2GRASP_TOR Summary of this function goes here
%   Detailed explanation goes here

%% Launch batch command
projectdirectory=GRASPTemplate.FolderPath;
work_directory=pwd;
cd(projectdirectory);
[status,cmdout]=system(['"' GRASPpath '" batch.gxp main.out main.log']);
% [status,cmdout]=system('grasp-analysis batch.gxp main.out main.log');
cd(work_directory);



end

