function [pattern] = readsphcutgrasp(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(filename);

freq=-1;

line=fgetl(fileID);%Comment
ph=1;
while(line~=-1)
    line=fgetl(fileID);%Header
    linenum=str2num(line);
    ICUT=linenum(6);
    if ICUT==1
        theta=linspace(linenum(1),linenum(1)+linenum(2)*(linenum(3)-1),linenum(3));
        pol=linenum(5);
        phi(ph)=linenum(4);       
        for th=1:length(theta)
            linenum=str2num(fgetl(fileID));
            E_1(th,ph)=linenum(1)+1i*linenum(2);
            E_2(th,ph)=linenum(3)+1i*linenum(4);
        end   
        ph=ph+1;
        line=fgetl(fileID);%Comment
    elseif ICUT==2
        phi=linspace(linenum(1),linenum(2)*(linenum(3)-1),linenum(3));
    end
end
field=zeros(length(theta),4,length(phi),length(freq));
if pol==3 || pol==2
    if ICUT==1
        field(:,1,:)=20*log10(abs(reshape(E_1,length(theta),1,length(phi))));
        field(:,2,:)=rad2deg(angle(reshape(E_1,length(theta),1,length(phi))));
        field(:,3,:)=20*log10(abs(reshape(E_2,length(theta),1,length(phi))));
        field(:,4,:)=rad2deg(angle(reshape(E_2,length(theta),1,length(phi))));
    elseif ICUT==2
        
    end
    pattern={field,theta',phi',freq,pol};
end
fclose(fileID);
end

