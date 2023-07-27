function GRASPTemplate = FFDualReflectorAFRMatlab2GRASP(GRASPTemplate,Array,Array_factor,GRASPpath,a,c)
%SINGLEREFLECTORMATLAB2GRASP Summary of this function goes here
%   Detailed explanation goes here

%% Matlab to GRASP
FFDualReflectorAFRMatlab2GRASP_TOR(GRASPTemplate,Array,Array_factor,a,c);
FFDualReflectorAFRMatlab2GRASP_TCI(GRASPTemplate,Array);
FFDualReflectorAFRMatlab2GRASP_Run(GRASPTemplate,GRASPpath);

resultsfilename='dual_cut.cut';
[E_co,E_cross] = FFDualReflectorAFRMatlab2GRASP_Read(GRASPTemplate,resultsfilename);

%% Save field
LPCP=GRASPTemplate.feed.LPCP;
pol=GRASPTemplate.feed.beta_0;

if LPCP == 1
    if pol == 0
        GRASPTemplate.FarfieldGRASP.E_H=E_co;
        GRASPTemplate.FarfieldGRASP.E_V=E_cross;
    elseif pol == pi/2
        GRASPTemplate.FarfieldGRASP.E_V=E_co;
        GRASPTemplate.FarfieldGRASP.E_H=E_cross;
    end
elseif LPCP == 2
    if pol == -1
        GRASPTemplate.FarfieldGRASP.E_LHCP=E_co;
        GRASPTemplate.FarfieldGRASP.E_RHCP=E_cross;
    elseif pol == 1
        GRASPTemplate.FarfieldGRASP.E_RHCP=E_co;
        GRASPTemplate.FarfieldGRASP.E_LHCP=E_cross;
    end
end

end

