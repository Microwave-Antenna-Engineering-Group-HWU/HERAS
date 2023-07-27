classdef SingleReflectorGRASPTemplate<GRASPtemplate
    %SINGLEREFLECTORGRASPTEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        reflector
        feed
        cut
        FarfieldGRASP
    end
    
    methods
        function obj = SingleReflectorGRASPTemplate(reflector,feed,cut,GRASPmethod)
            %SINGLEREFLECTORGRASPTEMPLATE Construct an instance of this class
            %   Detailed explanation goes here
%             obj.FolderPath=[pwd '\tempGRASP\'];
%             s = what('SingleReflector_batch');
%             copyfile(s.path,obj.FolderPath);
            obj.reflector=reflector;
            obj.feed=feed;
            obj.cut=cut;
            obj.TORfilename='SingleReflector.tor';
            obj.TCIfilename='SingleReflector.tci';
            obj.GRASPmethod=GRASPmethod;
        end

        function obj = CalculateTemplate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj = SingleReflectorMatlab2GRASP(obj);
        end
    end
end

