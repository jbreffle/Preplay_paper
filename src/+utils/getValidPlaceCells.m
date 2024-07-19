function validPlaceCell = getValidPlaceCells(PFresultsStruct_paramPoint, modelParam)
% validPlaceCell = utils.getValidPlaceCells(PFresultsStruct(ithParam1, ithParam2, ithNet), parameters);
%
% This function...
%
% Input:
%
%
% Ouput:
%

E_indices = PFresultsStruct_paramPoint.results{1}.E_indices;

PFmat_allTraj = [PFresultsStruct_paramPoint.results{1}.linfields{:}];
PFpeaks_allTraj = max(PFmat_allTraj(E_indices,:), [], 2);
validPlaceCell = PFpeaks_allTraj>=modelParam.minPeakRate;

end
