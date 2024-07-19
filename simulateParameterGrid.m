% sim_grid.m
%
% Runs simulations of both PFs and preplay across a grid of parameters. It
% saves the spikes of each as PFSpikeTimes_grid and preplaySpikeTimes_grid
%	-Each parameter point is simulated by the function run_param_sim()
% It then runs Bayesian decoding using those spikes.
%	-Decoding is performed by the function run_grid_decode
%
% Output:
%   Parameter structs:  modelParam, simParam
%   RNG data:           rngStruct_grid
%   Preplay spikes:     preplaySpikeTimes_grid{ithParam}{ithNet, ithEnv, ithCell}
%   PF spikes:          PFSpikeTimes_grid{ithParamSet}{ithNet, ithEnv, ithTrial, ithCell}
%
% Decode output:
%   paramOverride: analysis parameters different from the parameter struct
%   Preplay results:    resultsStruct(ithParam1, jthParam2, ithNet).results{1}
%       -fields: linfields, E_indices, netSeed, variedParamVals
%   PF results:         PFresultsStruct(ithParam1, jthParam2, ithNet).results
%       -fiels: various statistics and replaytrajectory has the decode data


%% Set up

addpath(fullfile("src"))
Config = utils.setConfig;
simSavePath = Config.simDataPath;


%% Simulation options

simParam.useTestParams = false; % 1 to perform a small-scale test simulation
simParam.useRandSeed = false; % false, to make sim_grid rng identical
simParam.testParamOption = [];
simParam.autoRun_sim_grid_decode = true; % Run run_grid_decode() right after all simulations complete
simParam.run_analysis_within_sim = 1; % 0 to skip all, 1 to skip decoding, 2 to run all

simParam.saveFlag = true; % If true, save simulation results
simParam.useDefaultSavePath = true; % 1 to select save destination, 0 to save in results dir

% simulateNetwork.m specific options:
simParam.runPFSim = [];	% Only needed for simulateNetwork.m
simParam.runPreplaySim = []; % Only needed for simulateNetwork.m
simParam.runDecode = []; % Only needed for simulateNetwork.m
simParam.plotResults = []; % Only needed for simulateNetwork.m
simParam.saveGvar = [];	% Only needed for simulateNetwork.m
simParam.saveVm = []; % Only needed for simulateNetwork.m

% Simulation options that are passed to the modelParam struct
simParam.calcPFScore = false; % If true, calculate PF score
simParam.nNets = 10; % Number of networks to simulate at each parameter point
simParam.nTrials_preplay = 1; % Number of preplay trials to simulate per network
simParam.t_max_preplay = 120; % Trial duration (s) for each preplay trial
simParam.nTrials_PF = 5; % Number of PF trials to simulate (each run across the track)
simParam.t_max_PF = 2; % Duration (s) for each PF trial
simParam.envIDs = [1:4];

simParam.nParPool = 4;
simParam.dispFlag_decode = 0;


% Simulate longer and more networks
runLongSims = false;
if runLongSims
    simParam.t_max_preplay = 300;
    simParam.nNets = 20;
    disp('Running long simulation')
end

% For simulating PFs with varied speed
runLongPFTrials = false;
if runLongPFTrials
    simParam.t_max_PF = 2; % 2, 4 seconds
    simParam.t_max_preplay = 5; % Very brief Preplay sim
    simParam.autoRun_sim_grid_decode = false; % Skip decoding of the brief preplay sim
end

%% Initialize parameters

simParam.paramSet = '7_22-v2' ;	% 'default', '7_22', '10_18', 'ISregMin', 'PFtuned1', 'PFtuned2'
disp(['Using ', simParam.paramSet, ' parameters'])
modelParam = initialize_parameters(simParam.paramSet, simParam);

connProbSet = 'standard'; % standard, low, high
disp(['Parameter connProbSet: ', connProbSet])

switch connProbSet
    case 'standard'
        % No changes
    case 'low'
        modelParam.conn_prob = 0.04;
        modelParam.del_G_syn_E_E = modelParam.del_G_syn_E_E * 2;
        disp('Overriding conn_prob to 0.04')
        modelParam = set_depedent_parameters(modelParam);
    case 'high'
        modelParam.conn_prob = 0.12;
        modelParam.del_G_syn_E_E = modelParam.del_G_syn_E_E * (2/3);
        disp('Overriding conn_prob to 0.12')
        modelParam = set_depedent_parameters(modelParam);
    otherwise
        error('Unknown connProbSet option')
end

% TODO: use a switch-case block to select these options

% Optional, for alternative feedforward input type
%modelParam.spatialInputType = 'stepped';
%modelParam.inputNSteps = 5;
%modelParam.inputNSteps = 3;
%modelParam.inputNSteps = 1;

% Three linearly varying cues
%modelParam.spatialInputType = 'linearTernary'; % linear, linearTernary, stepped

% Control simulation: without cluster correlations
%modelParam.clusterCorrs = false;


%% Set Up Grid Search Parameters
% Parameters must be a field set by initialize_parameters

gridSet = 'fiducial'; % mnc_clusters, Wee_Wie, dSRA_Wie, mnc_clusters_diag, fiducial
disp(['Parameter grid: ', gridSet])

switch gridSet
    case 'fiducial'
        variedParam(1).name = 'mnc';    % 1st parameter to be varied
        variedParam(1).range = [1.25];  %

        variedParam(2).name = 'clusters';	% 2nd parameter to be varied
        variedParam(2).range = [15];   %

    case 'mnc_clusters' % Fiducial is 1.25 and 15
        variedParam(1).name = 'mnc';    % 1st parameter to be varied
        variedParam(1).range = [1:0.25:3.0];  %

        variedParam(2).name = 'clusters';	% 2nd parameter to be varied
        variedParam(2).range = [5:5:30];   %

    case 'mnc_clusters_diag' % Fiducial is 1.25 and 15
        diagValues = [1, 2, 5, 10, 20, 30];
        variedParam(1).name = 'mnc';    % 1st parameter to be varied
        variedParam(1).range = diagValues;  %

        variedParam(2).name = 'clusters';	% 2nd parameter to be varied
        variedParam(2).range = diagValues;   %

    case 'Wee_Wie' % Fiducial is 220e-12 and 400e-12
        nLinspaceGrid = 7; % Number of parameters to test (each)

        variedParam(1).name = 'del_G_syn_E_E'; % 1st parameter to be varied
        baseVal1 = 220e-12;
        variedParam(1).range =  linspace( 0.75*baseVal1, 1.25*baseVal1, nLinspaceGrid);

        variedParam(2).name = 'del_G_syn_I_E'; % 2nd parameter to be varied
        baseVal2 = 400e-12;
        variedParam(2).range =  linspace( 0.75*baseVal2, 1.25*baseVal2, nLinspaceGrid);

        modelParam.clusters = 15;
        modelParam.mnc = 1.25;
        disp('Overriding parameters to manuscript values')
        assert(isequal(simParam.paramSet, '7_22-v2'))
        assert(isequal(modelParam.del_G_syn_I_E, modelParam.del_G_syn_E_I)) % Keep Wei fixed while Wie is varied
        modelParam = set_depedent_parameters(modelParam);

    case 'dSRA_Wie' % Fiducial is 3e-12 and 400e-12
        % nLinspaceGrid = 7; % Number of parameters to test (each)

        variedParam(1).name = 'del_G_sra'; % 1st parameter to be varied
        baseVal1 = 3e-12;
        variedParam(1).range = baseVal1.*logspace(0, 2, 5);

        variedParam(2).name = 'del_G_syn_I_E'; % 2nd parameter to be varied
        baseVal2 = 400e-12;
        variedParam(2).range =  linspace( 0.5*baseVal2, 2*baseVal2, 7);

        modelParam.clusters = 15;
        modelParam.mnc = 1.25;
        disp('Overriding parameters to manuscript values')
        assert(isequal(simParam.paramSet, '7_22-v2'))
        assert(isequal(modelParam.del_G_syn_I_E, modelParam.del_G_syn_E_I)) % Keep Wei fixed while Wie is varied
        modelParam = set_depedent_parameters(modelParam);

    otherwise
        error('Unknown gridSet option')
end


%% Parameters and hyperparameters for testing this code
if simParam.useTestParams
    disp(['Using testing params, testSize ' num2str(simParam.testParamOption)])
    simParam.run_analysis_within_sim = 2;
    simParam.useRandSeed = false;

    simParam.nNets = 2;
    simParam.nTrials_preplay = 1;
    variedParam(1).name = 'mnc'; % a parameter to be varied
    variedParam(1).range = [1.25]; %linspace(1, 8, 29); % set of values to test param2 at
    variedParam(2).name = 'clusters'; % a parameter to be varied
    variedParam(2).range = [10]; % set of values to test param2 at
    modelParam.t_max_preplay = 10;
    modelParam.envIDs = [1, 2];
    modelParam.shuffleIterations = 10;
    simParam.nParPool = 1;

    if simParam.testParamOption == 2
        % Slightly bigger test
        modelParam.t_max_preplay = 30;
        variedParam(1).range = [1.25, 2.0]; %linspace(1, 8, 29); % set of values to test param2 at
        variedParam(2).range = [10, 20]; % set of values to test param2 at
    elseif simParam.testParamOption == 3
        % Bigger test
        simParam.nNets = 4;
        modelParam.shuffleIterations = 500;
        modelParam.t_max_preplay = 120;
        variedParam(1).range = [1.25, 2.0]; %linspace(1, 8, 29); % set of values to test param2 at
        variedParam(2).range = [10, 20, 30]; % set of values to test param2 at
    end
end


%% Set up for sim loop

% In case any above changes affect dependent parameters
modelParam = set_depedent_parameters(modelParam);

% Combine into one parameter vector to pass to parfor function
parameterSets_vec = combvec(variedParam(:).range);

% Exclude cases where mnc>clusters
if isequal(variedParam(1).name, 'mnc') && isequal(variedParam(2).name, 'clusters')
    parameterSets_vec = parameterSets_vec(:,~[parameterSets_vec(1,:)>parameterSets_vec(2,:)]);
end

% Exclude cases where mnc>clusters
if isequal(variedParam(1).name, 'mnc') && isequal(variedParam(2).name, 'clusters')
    fullConnBoundary  = @(mnc) (mnc.^2)/modelParam.conn_prob;
    mncInd = find(strcmp({variedParam.name}, 'mnc'));
    mncVec = parameterSets_vec(mncInd,:);
    nClustInd = find(strcmp({variedParam.name}, 'clusters'));
    nClustVec = parameterSets_vec(nClustInd,:);
    badInds = fullConnBoundary(mncVec)<nClustVec;
    parameterSets_vec = parameterSets_vec(:,~badInds);
end

[~, simParam.grid_git_hash_string] = system('git rev-parse HEAD');
simParam.grid_working_dir = pwd;

if simParam.autoRun_sim_grid_decode
    disp('Will automatically run run_grid_decode() after simulations')
end
if simParam.run_analysis_within_sim==1
    disp('Will run event analysis within each simulation but not decoding')
end
if simParam.run_analysis_within_sim==2
    disp('Will run event decoding and analysis within each simulation')
end


%% Run grid search with spikes and spike stats returned

% Create save directory based on start date and time
simParam.sim_grid_name = datestr(now,'yyyy-mm-ddTHH-MM');
disp(simParam.sim_grid_name)
simParam.save_path = fullfile(simSavePath, ['grid_', simParam.sim_grid_name]);
if simParam.saveFlag
    if ~simParam.useDefaultSavePath
        simParam.save_path = uigetdir(simParam.save_path);
    end
    if ~isfolder(simParam.save_path)
        mkdir(simParam.save_path);
        disp('Using default save path')
    end
end
disp(['Save path: ', simParam.save_path])

% Starts parallel pool if needed
if simParam.nParPool>1
    disp(['Starting parpool with ', num2str(simParam.nParPool), ' workers'])
    if isempty(gcp('nocreate')) % create parpool if there isn't one
        parpool(simParam.nParPool);
    end
    current_par_pool = gcp('nocreate');
    if current_par_pool.NumWorkers~=simParam.nParPool
        delete(gcp('nocreate'))
        parpool(simParam.nParPool); % Create parpool if current one has wrong NumWorkers
    end
else
    disp('Using single worker')
end

% Store simParam
simParam.variedParam = variedParam;
simParam.parameterSets_vec = parameterSets_vec;

% Needed for tracking parfor progress
D = parallel.pool.DataQueue;
num_sets = size(parameterSets_vec, 2);
gridDispProgress(1, num_sets);
afterEach(D, @gridDispProgress);

resultsMatLinear = zeros(4, size(parameterSets_vec, 2));
resultsStructLinear = cell(1, size(parameterSets_vec, 2));
PFresultsStructLinear = cell(1, size(parameterSets_vec, 2));
preplaySpikeTimes_grid = cell(1, size(parameterSets_vec, 2));
PFSpikeTimes_grid = cell(1, size(parameterSets_vec, 2));
rngStruct_grid = cell(1, size(parameterSets_vec, 2));

if ~simParam.useRandSeed
    simParam.outsideSeed = 42;
    if modelParam.conn_prob==0.04
        simParam.outsideSeed = simParam.outsideSeed + 1;
    elseif modelParam.conn_prob==0.12
        simParam.outsideSeed = simParam.outsideSeed + 2;
    end
    disp(['Using rng seed ',  num2str(simParam.outsideSeed)])
else
    simParam.outsideSeed = randi(10^3);
    disp('Using random rng seed')
end
tic
if simParam.nParPool>1
    parfor ithParamSet = 1:size(parameterSets_vec, 2)
        rng(simParam.outsideSeed*ithParamSet, 'twister')
        mainSeed = randi(2^16);
        [PFSpikeTimes_grid{ithParamSet}, preplaySpikeTimes_grid{ithParamSet}, rngStruct_grid{ithParamSet},...
            ~, resultsStructLinear{ithParamSet}, PFresultsStructLinear{ithParamSet}] = ...
            run_param_sim(modelParam, simParam, ithParamSet, mainSeed);
        send(D, 1);
    end
else
    for ithParamSet = 1:size(parameterSets_vec, 2)
        rng(simParam.outsideSeed*ithParamSet, 'twister')
        mainSeed = randi(2^16);
        [PFSpikeTimes_grid{ithParamSet}, preplaySpikeTimes_grid{ithParamSet}, rngStruct_grid{ithParamSet},...
            ~, resultsStructLinear{ithParamSet}, PFresultsStructLinear{ithParamSet}] = ...
            run_param_sim(modelParam, simParam, ithParamSet, mainSeed);
        send(D, 1);
    end
end
simParam.simRunTime = toc;
disp(['Final simulation runtime: ', num2str(simParam.simRunTime/60/60), ' Hours'])
disp(['Finished simulation at ', char(datetime)])


%% Convert from linear structure to (ithParam1, ithParam2, ithNet)
resultsStruct = struct;
PFresultsStruct = struct;
for i = 1:size(parameterSets_vec, 2)
    structIndices = {};
    for ithParam = 1:size(variedParam, 2)
        structIndices{ithParam} = find(parameterSets_vec(ithParam,i)==variedParam(ithParam).range);
    end
    for ithNet = 1:simParam.nNets
        for k = 1:simParam.nTrials_preplay
            resultsStruct(structIndices{:}, ithNet, k).results = resultsStructLinear{i}{ithNet}{k};
        end
        PFresultsStruct(structIndices{:}, ithNet).results = PFresultsStructLinear{i}{ithNet};
    end
end


%% Save results
if simParam.saveFlag
    clear D h
    % Only save these as for now, until
    if simParam.useTestParams
        save(fullfile(simParam.save_path, ['grid_', simParam.sim_grid_name, 'tmp_ALL.mat']), '-v7.3')
    end
    if simParam.run_analysis_within_sim>0
        save(fullfile(simParam.save_path, ['grid_', simParam.sim_grid_name, 'results.mat']), 'resultsStruct', 'PFresultsStruct')
    end
    % Save spikes and parameter structurs
    save(fullfile(simParam.save_path, ['grid_', simParam.sim_grid_name, 'parameters.mat']), 'modelParam', 'simParam', 'rngStruct_grid');
    % Save PF spikes
    varinfo1=whos('PFSpikeTimes_grid'); saveopt1='';
    if (varinfo1.bytes >= 2^31); saveopt1='-v7.3'; end
    save(fullfile(simParam.save_path, ['grid_', simParam.sim_grid_name, 'spikes_PF.mat']), 'PFSpikeTimes_grid', saveopt1)
    % Save preplay spikes
    varinfo2=whos('preplaySpikeTimes_grid'); saveopt2='';
    if (varinfo2.bytes >= 2^31); saveopt2='-v7.3'; end
    save(fullfile(simParam.save_path, ['grid_', simParam.sim_grid_name, 'spikes_preplay.mat']), 'preplaySpikeTimes_grid', saveopt2)
end


%% Automatically run sim_grid_decode
if simParam.autoRun_sim_grid_decode
    disp('----------------------------------------')
    disp('Automatically starting run_grid_decode()')

    % Overwrite analysis parameters in the parameter struct
    paramOverride = struct;
    paramOverride.dispFlag_decode = 0;
    paramOverride.shuffleIterations = 100;
    %paramOverride.minPeakRate = 3;
    %paramOverride.nParPool = 1;
    %nParPool = 1;

    dataPath = simParam.save_path;
    nParPool = simParam.nParPool;
    [PFresultsStruct_post, resultsStruct_post] = run_grid_decode(dataPath, paramOverride=paramOverride, nParPool=nParPool);
    disp(['Finished decode at ', char(datetime)])

else
    disp('Did not run simulation decoding')
    disp('Did not add results to gridSimNames.json')
end

