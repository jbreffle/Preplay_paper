function eventRippleIds = getValidEventIDs(PFresultsStruct_paramPoint, resultsStruct_paramPoint, modelParam)
% eventRippleIds = utils.getValidEventIDs(PFresultsStruct(ithParam1, ithParam2, ithNet), resultsStruct(ithParam1, ithParam2, ithNet), parameters);
%
% Criterion for preplay events are stricter than criterion for detecting
% PBE events.
%
% This function returns the indices of the PBE events that are valid
% preplay events.
%
% To go from ripple to event: must meet min duration, min cells
% participating, and participating cells must meet min peak PF rate.
% Peak PF rate is taken from all simulated trajectories.
%
% Input:
%   Results structures for particular net at particular parameter point
%   model parameters
%
% Ouput:
%   Index of the PBE events that are valid preplay events

validPlaceCell = utils.getValidPlaceCells(PFresultsStruct_paramPoint, modelParam);

rippleEventLengths = resultsStruct_paramPoint.results.ripple.endtime-resultsStruct_paramPoint.results.ripple.starttime;
validEventLength = rippleEventLengths>(modelParam.minEventDur/1000);

%rippleCellCounts = round(resultsStruct(ithParam1, ithParam2, ithNet).results.eventParticipation*parameters.n_E)';
%rippleCellCounts = sum(~isnan(resultsStruct(ithParam1, ithParam2, ithNet).results.ranksVec))';
rippleCellCounts = sum(~isnan(resultsStruct_paramPoint.results.ranksVec) & validPlaceCell)';
validEventCellCount = rippleCellCounts>=modelParam.cellcountthresh;

eventRippleIds = find(validEventLength&validEventCellCount);

end
