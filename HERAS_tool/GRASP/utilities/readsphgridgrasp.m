function [E_1,E_2,u,v] = readsphgridgrasp(filename)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


%% Read header
fileID = fopen(filename);

line=fgetl(fileID);
while(~strcmp(line,'++++'))
    line=fgetl(fileID);
end

for i=1:3
    fgetl(fileID);
end

line=fgetl(fileID);%Comment
linenum=str2num(line);

u_s=linenum(1);
v_s=linenum(2);
u_e=linenum(3);
v_e=linenum(4);

line=fgetl(fileID);%Comment
linenum=str2num(line);

u_single=linspace(u_s,u_e,linenum(1));
v_single=linspace(v_s,v_e,linenum(2));

[u,v]=meshgrid(u_single,v_single);

%% Read data
data=textscan(fileID,'%f %f %f %f');

nan_index=find(isnan(data{4})==1);
N=sum(nan_index);

if N==0
    E_1=reshape(data{1,1}(:)+1j*data{1,2}(:),size(u')).';
    E_2=reshape(data{1,3}(:)+1j*data{1,4}(:),size(u')).';
else
    E_1=nan(linenum(2),linenum(1));
    E_2=E_1;
    for nv=1:linenum(2)
        column_start=data{1,1}(nan_index(nv));
        N_columns=data{1,2}(nan_index(nv));
        E_1(nv,column_start+(0:N_columns-1))=reshape(data{1,1}(nan_index(nv)+(1:N_columns))+1j*data{1,2}(nan_index(nv)+(1:N_columns)),[1,N_columns]);
        E_2(nv,column_start+(0:N_columns-1))=reshape(data{1,3}(nan_index(nv)+(1:N_columns))+1j*data{1,4}(nan_index(nv)+(1:N_columns)),[1,N_columns]);
    end
end

fclose(fileID);

end

