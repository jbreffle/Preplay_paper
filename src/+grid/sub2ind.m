function ind = sub2ind(ithParam1, ithParam2, simParam)
% ind = grid.sub2ind(ithParam1, ithParam2, simParam);
%
% Convert parameter grid linear index to 2D subscripts
% Convert linear ind to 2d Ind
param1Val = simParam.variedParam(1).range(ithParam1);
param2Val = simParam.variedParam(2).range(ithParam2);
ind = find(ismember( ...
    simParam.parameterSets_vec', [param1Val; param2Val]', 'rows' ...
    ));
%assert(~isempty(exampleLinearParamInd), "grid.ind2sub() error")
end
