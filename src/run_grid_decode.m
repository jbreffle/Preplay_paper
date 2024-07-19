function [PFresultsStruct, resultsStruct] = run_grid_decode(data_path, varargin)
% run_grid_decode(simParam.save_path)
%
% Code to perform bayesian decoding on the spike-time cell arrays produced
% by simulateParameterGrid.m. Adapted from the script run_grid_decode
%
% Inputs:
%   - data_path: a string of the location where simulation results were saved
%
% Outputs:
%   - PFresultsStruct: PF results struct
%   - resultsStruct: preplay results struct
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj,   'data_path',                    @(x) isstring(x)||ischar(x) )
addParameter(inputObj,	'plotInterimResults',   false,	@islogical);
addParameter(inputObj,	'plotSummaryResults',	true,	@islogical);
addParameter(inputObj,	'nParPool',             1,      @isnumeric);
addParameter(inputObj,	'paramOverride',        1,      @isstruct);
parse(inputObj, data_path, varargin{:});
p = inputObj.Results;
data_path = char(data_path);


%% Set up

Config = utils.setConfig;
% TODO: switch to using Config struct

[~, sim_grid_folder_name] = fileparts(data_path);
decode_name = datestr(now,'yyyy-mm-ddTHH-MM');

% Load data
load(fullfile(data_path, [sim_grid_folder_name, 'parameters.mat']), 'simParam') % Load the necessary parameter structures
load(fullfile(data_path, [sim_grid_folder_name, 'spikes_PF.mat']), 'PFSpikeTimes_grid'); % Load all spikes now, pass subsets to run_param_decode()
load(fullfile(data_path, [sim_grid_folder_name, 'spikes_preplay.mat']),  'preplaySpikeTimes_grid'); % Load all spikes now, pass subsets to run_param_decode()

disp(['data_path ', data_path])
disp(['sim_grid_folder_name ', sim_grid_folder_name])
disp(['decode_name ', decode_name])

disp('Overriding these parameters:')
disp(fieldnames(p.paramOverride)')


%% Start parpool pool if needed
if p.nParPool>1
    disp(['Using parpool with ', num2str(p.nParPool), ' workers'])
    if isempty(gcp('nocreate')) % create parpool if there isn't one
        parpool(p.nParPool);
    end
    current_par_pool = gcp('nocreate');
    if current_par_pool.NumWorkers~=p.nParPool
        delete(gcp('nocreate'))
        parpool(p.nParPool); % Create parpool if current one has wrong NumWorkers
    end
else
    delete(gcp('nocreate'))
    disp('Using single worker')
end

% Needed for tracking parfor progress
D = parallel.pool.DataQueue;
num_sets = size(simParam.parameterSets_vec, 2);
gridDispProgress(1, num_sets);
afterEach(D, @gridDispProgress);


%% Perform analysis, looping over parameter sets
%paramSetInds = combvec([1:numel(simParam.variedParam(1).range)], [1:numel(simParam.variedParam(2).range)])';
[param1IndVec, ~] = find([simParam.variedParam(1).range==simParam.parameterSets_vec(1,:)']');
[param2IndVec, ~] = find([simParam.variedParam(2).range==simParam.parameterSets_vec(2,:)']');
paramSetInds = [param1IndVec, param2IndVec];

PFresultsStruct_linear = struct;
resultsStruct_linear = struct;

tic
if ~isempty(gcp('nocreate'))
    parfor ithParamSet = 1:size(simParam.parameterSets_vec, 2)
        [PFresults, allResults] = run_param_decode(PFSpikeTimes_grid{ithParamSet}, preplaySpikeTimes_grid{ithParamSet}, data_path, paramSetInds, ithParamSet, ...
            plotInterimResults=p.plotInterimResults, paramOverride=p.paramOverride);
        PFresultsStruct_linear(ithParamSet).net = PFresults;
        resultsStruct_linear(ithParamSet).net = allResults;
        send(D, 1);
    end
else
    for ithParamSet = 1:size(simParam.parameterSets_vec, 2)
        [PFresults, allResults] = run_param_decode(PFSpikeTimes_grid{ithParamSet}, preplaySpikeTimes_grid{ithParamSet}, data_path, paramSetInds, ithParamSet, ...
            plotInterimResults=p.plotInterimResults, paramOverride=p.paramOverride);
        PFresultsStruct_linear(ithParamSet).net = PFresults;
        resultsStruct_linear(ithParamSet).net = allResults;
        send(D, 1);
    end
end
decodeRunTime_seconds = toc;
disp(['Decode RunTime: ', num2str(decodeRunTime_seconds/60/60), ' Hours'])


%% Un-linearize the structs
PFresultsStruct = struct;
resultsStruct = struct;
for ithParamSet = 1:size(simParam.parameterSets_vec, 2)
    ithParam1 = paramSetInds(ithParamSet,1);
    ithParam2 = paramSetInds(ithParamSet,2);
    for ithNet = 1:simParam.nNets
        PFresultsStruct(ithParam1, ithParam2, ithNet).results{1} = PFresultsStruct_linear(ithParamSet).net{ithNet}{1};
        for ithTest = 1:simParam.nTrials_preplay
            resultsStruct(ithParam1, ithParam2, ithNet, ithTest).results = resultsStruct_linear(ithParamSet).net{ithNet}{ithTest};
        end
    end
end


%% Save results, consistent with previous simulation data structures
varinfo1=whos('resultsStruct'); varinfo2=whos('PFresultsStruct');
saveopt=''; % use default of save()
if (varinfo1.bytes >= 2^31) || (varinfo2.bytes >= 2^31)
    saveopt='-v7.3'; % use version that allows large files
end
if simParam.saveFlag
    paramOverride = p.paramOverride;
    fullDecodeName = ['grid_', simParam.sim_grid_name, 'decode_', decode_name];
    save( fullfile(data_path, fullDecodeName) , ...
        'resultsStruct', 'PFresultsStruct', 'paramOverride', 'decodeRunTime_seconds', saveopt);

    % Write decode name to gridSimNames.json
    utils.writeToGridSimNameJson(fullDecodeName)

else
    disp('    NOT SAVING RESULTS    ')
end


%% Plot place fields and preplay summary stats
if p.plotSummaryResults
    % Place field...
    % Summary stats ...
end


%% Plot preplay decoding results
if p.plotSummaryResults
    ithEnv = 1; ithNet=1; ithTest=1;
    ithParam1 = 1; ithParam2 = 1;
    replaytrajectory = resultsStruct(ithParam1, ithParam2, ithNet, ithTest).results.replaytrajectory;

    if ~isempty(replaytrajectory) && ~isempty(replaytrajectory.pvalue)
        pvals = replaytrajectory.pvalue(:,ithEnv);
        figure; histogram(pvals, 10)

        allshuff_rvals = vertcat(replaytrajectory.shuffle_rsquare{:});
        allshuff_rvals = allshuff_rvals(:,1); % take just forward traj
        rvals_preplay = replaytrajectory.rsquare(:,1);
        rvals_shuffle = vertcat(allshuff_rvals{:,1});
        [~,P,KSSTAT] = kstest2(rvals_preplay, rvals_shuffle);
        disp(['Decode pval ', num2str(P)])

        figure; hold on;
        ecdf(rvals_preplay)
        ecdf(rvals_shuffle)
        legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
    end
end


end % end of run_grid_decode()