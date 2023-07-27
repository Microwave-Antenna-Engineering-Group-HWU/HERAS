classdef GRASPtemplate<matlab.mixin.Copyable
    %GRASPTEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FolderPath
        TORfilename
        TCIfilename
        OUTPUTfilename
        GRASPmethod(1,:) char {mustBeMember(GRASPmethod,{'po','ptd','po_plus_ptd'})} = 'po'
    end
    
    methods
        function obj = GRASPtemplate(~)
            %GRASPTEMPLATE Construct an instance of this class
            %   Detailed explanation goes here
        end
        
    end
end

