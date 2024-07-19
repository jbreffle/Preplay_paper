%% This script performs repeated-decodings of a given simulation with varied decode parameters
% varyDecodeParam.m
%
% run_grid_decode.m performs decoding across a parameter grid.
%   It calls run_param_decode to decode each parameter, and this calls
%       run_net_decode.m to decode each network.
%
% Outline:
% - Load a previously completed simulation
% - Select a particular parameter point to analyze
% - Loop over values of a decode parameter (such as fraction of cells)
%   - Call run_param_decode
%   - Only save subset of its output? The replaytrajectory field might be
%   all we care about (but the other fields are small, so might be no harm
%   in keeping them)
%   - Need to call run_param_decode multiple times with different rng to
%   get a distribution of results for each decode parameter
%
%
% Decode data structure produced by run_grid_decode:
% ithParam1=1; ithParam2=1; ithNet=1;
% resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory
%
% The function run_param_decode.m runs analysis of one parameter point from
% results produced by simulateParameterGrid.m
%
% The function run_grid_decode.m performs the analysis for an entire grid
% by calling run_param_decode.m for each paramter in the grid.
%
% See gridAnalysisSpiking.m for an example of loading and using
% parameter-specific spike data
%
% See simulateParameterGrid.m approach to allowing arbitrary parameters to be varied in
% the decoding parameter value loop
%
% Need to over write the parameters
%    modelParam.downSampleCells = false; % true/false flag to use downsampling
%    modelParam.downSampledFrac = 0.1;  % fraction of cells to use in decode, if downSampleCells==true


addpath(fullfile("..", "src"))
Config = utils.setConfig;
% TODO: switch to using Config struct

% Automatically get data location
userDir = char(java.lang.System.getProperty('user.home'));
if ispc; boxDir = ['Box'];
elseif ismac; boxDir = ['Library' filesep 'CloudStorage' filesep 'Box-Box'];
else; disp('error');
end
newSimDir = ['Data' filesep 'RandNet-Data' filesep 'new-structure'];
dataPath = [userDir filesep boxDir filesep newSimDir filesep];
oldDataPath = [userDir filesep boxDir filesep 'Data' filesep 'RandNet-Data' filesep 'RandNet param sweeps'];


%% Load simulation data

decodeFileID = "mainPaperGrid"; % Select which decoding file result ID to analyze

% 216 is main result for manuscript
% 224 is similar, but will run more quickly due to single environment for testing this code
analysisParam.simName = grid.getDecodeFilename(decodeFileID);

tmp = split(analysisParam.simName, 'decode');
simName = tmp{1};
load([dataPath filesep simName filesep simName 'parameters.mat'])
load([dataPath filesep simName filesep analysisParam.simName '.mat'])
load([dataPath filesep simName filesep simName 'spikes_preplay.mat'])
load([dataPath filesep simName filesep simName 'spikes_PF.mat'])

analysisParam.ithSimParam1 = 2;
analysisParam.ithSimParam2 = 3;

param1Val = simParam.variedParam(1).range(analysisParam.ithSimParam1);
disp(['Param1=', num2str(analysisParam.ithSimParam1), ': ', simParam.variedParam(1).name, '=', num2str(param1Val)])
param2Val = simParam.variedParam(2).range(analysisParam.ithSimParam2);
disp(['Param2=', num2str(analysisParam.ithSimParam2), ': ', simParam.variedParam(2).name, '=', num2str(param2Val)])

linearParamInd = find(all([param1Val; param2Val] == simParam.parameterSets_vec, 1));


PFSpikeTimes_paramPoint = PFSpikeTimes_grid{linearParamInd};
preplaySpikeTimes_paramPoint = preplaySpikeTimes_grid{linearParamInd};
clear preplaySpikeTimes_grid PFSpikeTimes_grid

%% Analysis parameters

%analysisParam.nReplicates = 6;	% Number of times to repeate each decoding analysis
analysisParam.nReplicates = 25;	% Number of times to repeate each decoding analysis
analysisParam.nWorkers = 6;     % Number of parfor workers

analysisParam.saveResults = true;

% Note: 375*0.03125 = 11.72
% Note: 10^-1.5 = 0.0316
analysisParam.variedParam.name = 'downSampledFrac';
%analysisParam.variedParam.range = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125];
analysisParam.variedParam.range = logspace(0, -1.5, 13);
%analysisParam.variedParam.range = [1.0, 0.25];
%analysisParam.variedParam.range = [1.0, 0.5, 0.25];
modelParam.downSampleCells = true;

paramOverride.shuffleIterations = 100;
% paramOverride.shuffleIterations = 10
paramOverride.downSampleCells = true;


%% Starts parallel pool if needed
if analysisParam.nWorkers>1
    disp(['Starting parpool with ', num2str(analysisParam.nWorkers), ' workers'])
    if isempty(gcp('nocreate')) % create parpool if there isn't one
        parpool(analysisParam.nWorkers);
    end
    current_par_pool = gcp('nocreate');
    if current_par_pool.NumWorkers~=analysisParam.nWorkers
        delete(gcp('nocreate'))
        parpool(analysisParam.nWorkers); % Create parpool if current one has wrong NumWorkers
    end
else
    disp('Using single worker')
end

[~, analysisParam.git_hash_string] = system('git rev-parse HEAD');
analysisParam.running_file = mfilename('fullpath');

%% Perform decoding, looping over decode parameter values

[param1IndVec, ~] = find([simParam.variedParam(1).range==simParam.parameterSets_vec(1,:)']');
[param2IndVec, ~] = find([simParam.variedParam(2).range==simParam.parameterSets_vec(2,:)']');
allParamSetInds = [param1IndVec, param2IndVec];
assert( isequal(allParamSetInds(linearParamInd,:), [analysisParam.ithSimParam1, analysisParam.ithSimParam2]) )

decode_data_path = [dataPath simName];
%decodeParamSetInds = [analysisParam.ithSimParam1, analysisParam.ithSimParam2];
%decodeLinearParamInd = 1;

% Initialize the main output structures of the analysis
PFresultsStruct_variedDecode = struct;
resultsStruct_variedDecode = struct;

% modelParam.downSampleCells must be true for modelParam.downSampledFrac to affect the analysis
if isequal(analysisParam.variedParam.name,  'downSampledFrac')
    assert(paramOverride.downSampleCells)
end

tic
for ithParam = 1:numel(analysisParam.variedParam.range)

    % Update the parameter that is looped over
    paramOverride.(analysisParam.variedParam.name) = analysisParam.variedParam.range(ithParam);

    % Loop over random replicates of the analysis
    parfor ithReplicate = 1:analysisParam.nReplicates
        disp(['Starting replicate ', num2str(ithReplicate), ' of ithParam ', num2str(ithParam)])

        cantorPairingSeed = 0.5*(ithParam+ithReplicate)*(ithParam+ithReplicate+1)+ithReplicate;
        %cantorPairingSeed = 27330*ithReplicate;

        [PFresults, allResults] = run_param_decode(PFSpikeTimes_paramPoint, preplaySpikeTimes_paramPoint, decode_data_path, allParamSetInds, linearParamInd, ...
            plotInterimResults=false, paramOverride=paramOverride, mainSeed=cantorPairingSeed);
        PFresultsStruct_variedDecode(ithParam, ithReplicate).net = PFresults;
        resultsStruct_variedDecode(ithParam, ithReplicate).net = allResults;
    end

    disp(['Finished ithParam ', num2str(ithParam), '. Runtime: ', num2str(toc/60), 'min'])
end
runTime = toc;
disp(['Runtime: ', num2str(runTime/60), 'min'])


%% Save results

if analysisParam.saveResults

    % Remove un-necessary fields that are too large
    resultsStruct_variedDecode = utils.removeExtraDecodeStructFields(resultsStruct_variedDecode);

    % Save results
    fileName = [simName, 'variedDecode', datestr(now,'yyyy-mm-ddTHH-MM')];
    save_path = fullfile(dataPath, simName, fileName);
    save(save_path, 'PFresultsStruct_variedDecode', 'resultsStruct_variedDecode', 'analysisParam', 'modelParam', 'runTime', 'paramOverride')

    disp(['Saved file ', fileName])
else
    disp('DID NOT SAVE RESULTS')
end

