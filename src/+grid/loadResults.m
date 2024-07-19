function [modelParam, simParam, resultsStruct, PFresultsStruct, varargout] = loadResults(decodeName)
% [modelParam, simParam, resultsStruct, PFresultsStruct, preplayGridSpikes, pfGridspikes] = grid.loadResults(decodeName)
%
% Returns preplay spikes as optional 5th output
% Returns place field spikes as optional 6th output


%% Set up

Config = utils.setConfig;
dataPath = Config.simDataPath;


%% Load files

splitFilename = split(decodeName, 'decode');
simName = splitFilename{1};

% modelParam, rngStruct_grid, simParam
loadedParams = load(fullfile(dataPath, simName, [simName 'parameters.mat']));

% PFresultsStruct, resultsStruct, decodeRunTime_seconds, paramOverride
loadedDecode = load(fullfile(dataPath, simName, [decodeName '.mat']));

modelParam = loadedParams.modelParam;
rngStruct_grid = loadedParams.rngStruct_grid;
simParam = loadedParams.simParam;

PFresultsStruct = loadedDecode.PFresultsStruct;
resultsStruct = loadedDecode.resultsStruct;
decodeRunTime_seconds = loadedDecode.decodeRunTime_seconds;
paramOverride = loadedDecode.paramOverride;

% load and return preplaySpikeTimes_grid and PFSpikeTimes_grid if requested
% First preplay, then PF
if nargout > 4
    % preplaySpikeTimes_grid
    loadedPreplaySpikes = load(fullfile(dataPath, simName, [simName 'spikes_preplay.mat']));
    preplaySpikeTimes_grid = loadedPreplaySpikes.preplaySpikeTimes_grid;
    varargout{1} = preplaySpikeTimes_grid;
end
if nargout == 6
    % PFSpikeTimes_grid
    loadedPfSpikes = load(fullfile(dataPath, simName, [simName 'spikes_PF.mat']));
    PFSpikeTimes_grid = loadedPfSpikes.PFSpikeTimes_grid;
    varargout{2} = PFSpikeTimes_grid;
end

end
