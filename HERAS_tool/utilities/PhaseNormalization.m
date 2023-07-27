function [PhaseNormalizedMatrix] = PhaseNormalization(InputMatrix,TargetPhase,TargetPhaseIndex)
% This function recives an NxM phase matrix and normalizes it so that
% PhaseMatrix(TargetPhaseIndex(i),i)=TargetPhase(i). 

[N,M]=size(InputMatrix);
N_TargetPhase=length(TargetPhase);
N_TargetPhaseIndex=length(TargetPhaseIndex);
if (N_TargetPhase~=N_TargetPhaseIndex)||(N_TargetPhase~=M)
    error('Wrong input size');
end

PhaseNormalizedMatrix=zeros(N,M);

parfor m=1:M
    PhaseNormalizedMatrix(:,m)=InputMatrix(:,m).*exp(1i*(TargetPhase(m)-angle(InputMatrix(TargetPhaseIndex(m),m))))
end
