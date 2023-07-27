classdef DualReflectorGRASPTemplate<GRASPtemplate
    %SINGLEREFLECTORGRASPTEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        reflector
        subreflector
        feed
        cut
        FarfieldGRASP
    end
    
    methods
        function obj = DualReflectorGRASPTemplate(reflector,subreflector,feed,cut,GRASPmethod)
            %SINGLEREFLECTORGRASPTEMPLATE Construct an instance of this class
            %   Detailed explanation goes here
%             obj.FolderPath=[pwd '\tempGRASP\'];
%             s = what('SingleReflector_batch');
%             copyfile(s.path,obj.FolderPath);
            obj.reflector=reflector;
            obj.subreflector=subreflector;
            obj.feed=feed;
            obj.cut=cut;
            obj.TORfilename='DualReflector_batch.tor';
            obj.TCIfilename='DualReflector_batch.tci';
            obj.GRASPmethod=GRASPmethod;
        end

        function obj = CalculateTemplate(obj,GRASPpath,a,c)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj = FFDualReflectorAFRMatlab2GRASP(obj,obj.feed,1,GRASPpath,a,c);
        end
    end
end

