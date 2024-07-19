function [ithParam1, ithParam2] = ind2sub(linearParamInd, simParam)
% [ithParam1, ithParam2] = grid.ind2sub(linearParamInd, simParam);
%
% Convert parameter grid linear index to 2D subscripts
ithParam1Val = simParam.parameterSets_vec(1,linearParamInd);
ithParam2Val = simParam.parameterSets_vec(2,linearParamInd);
ithParam1 = find(simParam.variedParam(1).range==ithParam1Val);
ithParam2 = find(simParam.variedParam(2).range==ithParam2Val);
%assert(~isempty([ithParam1, ithParam2]), "grid.ind2sub() error")
end
