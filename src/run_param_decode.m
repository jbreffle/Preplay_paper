function [PFresults, allResults] = run_param_decode(PFSpikeTimes_paramPoint, preplaySpikeTimes_paramPoint, data_path, paramSetInds, ithParamSet, varargin)
% [PFresults, allResults] = run_param_decode(PFSpikeTimes_paramPoint, preplaySpikeTimes_paramPoint, data_path, paramSetInds, ithParamSet)
%
% Runs analysis and Bayesian decoding for all networks at one parameter
% point from a parameter grid simulation.
%
% Inputs:
%   - PFSpikeTimes_paramPoint: all PF spikes for one parameter point
%   - preplaySpikeTimes_paramPoint: all preplay spikes for one parameter point
%   - data_path: The location that the data is saved in
%   - paramSetInds: The vector of parameter grid parameter changes
%   - ithParamSet: The index of paramSetInds that corresponds to the
%   current call to run_param_decode()
%
% Outputs:
%   - PFresults: the Place Field analysis struct
%   - allResults: the preplace analysis struct
%
% Optional:
%   - plotInterimResults, paramOverride
%
% Loads: modelParam, simParam, and rngStruct_grid from the parameters.mat file

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj,   'PFSpikeTimes_paramPoint',      @iscell)
addRequired(inputObj,   'preplaySpikeTimes_paramPoint',	@iscell)
addRequired(inputObj,   'data_path',            @ischar)
addRequired(inputObj,   'paramSetInds',         @isnumeric)
addRequired(inputObj,   'ithParamSet',          @isnumeric)
addParameter(inputObj,   'plotInterimResults',  false,  @islogical)
addParameter(inputObj,   'paramOverride',       [],     @isstruct)
addParameter(inputObj,   'mainSeed',      [],     @isnumeric)
parse(inputObj, PFSpikeTimes_paramPoint, preplaySpikeTimes_paramPoint, data_path, paramSetInds, ithParamSet, varargin{:});
p = inputObj.Results;

[~, sim_grid_folder_name] = fileparts(data_path);


%% Load params and set-up for analysis
load(fullfile(data_path, [sim_grid_folder_name, 'parameters.mat']), 'modelParam', 'simParam', 'rngStruct_grid') % Load the 3 parameter structures

% Save all results in these fields
day = 1; eprun = 2; tet = 1;

% Decoding function options
decodeEpoch=2;
savedata_withinDecoding = false;
figopt_withinDecoding = 0;

% Overwrite analysis parameters in the parameter struct
if ~isempty(p.paramOverride)
    for field_name = fieldnames(p.paramOverride)'
        modelParam.(field_name{1}) = p.paramOverride.(field_name{1});
    end
end


%% Run analyses
% ithParam1 = paramSetInds(ithParamSet,1);
% ithParam2 = paramSetInds(ithParamSet,2);

PFresults = cell(1, simParam.nNets); % Initialize net's struct
allResults = cell(1, simParam.nNets); % Initialize net's struct

for ithNet = 1:simParam.nNets

    %% Recreate network
    E_indices = rngStruct_grid{ithParamSet}.net(ithNet).E_indices;
    netSeed = rngStruct_grid{ithParamSet}.net(ithNet).netSeed;
    netParams=modelParam;
    netParams.envIDs = modelParam.envIDs;% pfsim.envIDs;
    for ithParam = 1:size(simParam.variedParam, 2)
        netParams.(simParam.variedParam(ithParam).name) = simParam.variedParam(ithParam).range(paramSetInds(ithParamSet,ithParam));
    end
    netParams = set_depedent_parameters(netParams);
    network = create_network(netParams, 'seed', netSeed);
    if ~all(network.E_indices==E_indices); error('Incorrect network'); end


    %% Place Field analysis
    % Extract all Place Field spikes of this network
    netSpikeTimes = reshape(PFSpikeTimes_paramPoint(ithNet,:,:,:), [modelParam.nEnvironments, modelParam.nTrials_PF, modelParam.n]);
    %opS = false(pfsim.nEnvironments, pfsim.nTrials, parameters.n, pfsim.t_steps);
    opS = false(modelParam.n, numel(modelParam.t_PF), modelParam.nEnvironments, modelParam.nTrials_PF);
    for ithEnv = 1:modelParam.nEnvironments
        for ithTrial = 1:modelParam.nTrials_PF
            for ithCell = 1:modelParam.n %E_indices
                % cellSpikeInds = round(netSpikeTimes{ithCell}./parameters.dt);
                cellSpikeInds = round(netSpikeTimes{ithEnv, ithTrial, ithCell}./modelParam.dt);
                opS(ithCell, cellSpikeInds, ithEnv, ithTrial) = 1;
            end
        end
    end
    if p.plotInterimResults % Plot mean Place field sim rates
        figure; plot(movmean( squeeze(mean(opS(:,:,ithEnv,:), [1,4])), 100))
    end

    % Recreate linfields
    linfields = calculate_linfields(opS, modelParam);
    [PFmat, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);

    % Store PF results in struct
    PFresults{ithNet}{1}.linfields = PFmat;
    PFresults{ithNet}{1}.E_indices = network.E_indices;
    PFresults{ithNet}{1}.netSeed = netSeed;
    PFresults{ithNet}{1}.variedParamVals = simParam.parameterSets_vec(:,ithParamSet);


    %% Preplay analysis
    for ithTest = 1:simParam.nTrials_preplay

        % Get preplay sim spike times
        netSpikeTimes = squeeze(preplaySpikeTimes_paramPoint(ithNet, ithTest, :));
        EspikeMat = false(numel(E_indices), modelParam.t_steps_preplay);
        for ithECell = 1:numel(E_indices) %E_indices
            ithECellInd = E_indices(ithECell);
            cellSpikeInds = round(netSpikeTimes{ithECellInd}./modelParam.dt);
            EspikeMat(ithECell,cellSpikeInds) = 1;
        end
        if p.plotInterimResults % Plot mean Place field sim rates
            figure; plot(movmean( mean(EspikeMat, 1), 100))
        end
        % Identify PBEs
        [trialResults] = detect_PBE(EspikeMat, modelParam);

        % Accumulate trialResults over ithTests
        if ithTest == 1 % append trialResults struct to network results struct
            network_spike_sequences = trialResults;
        else
            network_spike_sequences = [network_spike_sequences, trialResults];
        end
        % Analyze and save summary statistics
        % Calculate mean n spikes within events
        overallmeanspikes = 0;
        for ithEvent = 1:size(trialResults.events, 1)
            eventSpikes = EspikeMat(:,trialResults.events(ithEvent,1):trialResults.events(ithEvent,2));
            nSpikes = sum(eventSpikes, 2);
            meannSpikes = sum(nSpikes)/sum(nSpikes>0);
            overallmeanspikes = overallmeanspikes+ meannSpikes;
        end
        overallmeanspikes = overallmeanspikes./size(trialResults.events, 1);
        % Overall simulation statistics
        allResults{ithNet}{ithTest}.ithInit = ithTest;
        allResults{ithNet}{ithTest}.numEvents = numel(network_spike_sequences(ithTest).event_lengths); % number of detected events
        allResults{ithNet}{ithTest}.fracFire =  mean(sum(EspikeMat, 2)>0); % Fraction of cells that fire at all during simulation
        allResults{ithNet}{ithTest}.frac_participation = mean([network_spike_sequences(ithTest).frac_spike{:}]); % mean fraction of cells firing per event
        allResults{ithNet}{ithTest}.meanRate = mean(sum(EspikeMat, 2)/ modelParam.t_max_preplay); % mean over cells' average firing rate
        allResults{ithNet}{ithTest}.stdRate = std(sum(EspikeMat, 2)/ modelParam.t_max_preplay); % STD over cells' average firing rate
        allResults{ithNet}{ithTest}.meanCellnEventSpikes = overallmeanspikes; % mean n event spikes over cells that spike
        % Stats for each detected event
        allResults{ithNet}{ithTest}.eventLength = network_spike_sequences(ithTest).event_lengths; % duration in seconds of all detected events
        allResults{ithNet}{ithTest}.eventParticipation = [network_spike_sequences(ithTest).frac_spike{:}]; % fraction of cells that fired in each event
        % Save ranks_vec
        allResults{ithNet}{ithTest}.ranksVec = network_spike_sequences(ithTest).ranks_vec;


        %% Bayesian Decoding, for each preplay trial

        allTrajEPFs = [];
        for iTraj = 1:numel(PFmat)
            allTrajEPFs = [allTrajEPFs, PFmat{iTraj}(network.E_indices,:)];
        end
        arePlaceCells = max(allTrajEPFs, [], 2)>modelParam.minPeakRate;
        %disp(['There are ', num2str(sum(arePlaceCells)), ' place cells'])

        % Calculate and save preplay Bayesian decoding
        % Only decode if there are events and place cells
        if ~isempty(trialResults.events) && any(arePlaceCells)
            if ~isempty(p.mainSeed) % Over-ride mainSeed if passed in as optional parameter
                mainSeed = p.mainSeed;
            else
                mainSeed = rngStruct_grid{ithParamSet}.mainSeed;
            end
            if simParam.run_analysis_within_sim>1
                decodeFuncSeed = rngStruct_grid{ithParamSet}.net(ithNet).decodeFunc(ithTest);
            else
                decodeFuncSeed = ithTest*ithNet*1000*mainSeed;
            end
            % Run decode
            replaytrajectory = run_net_decode(modelParam, network, trialResults, linfields, EspikeMat, day=day, ep=decodeEpoch, figOpt=figopt_withinDecoding, decodeFuncSeed=decodeFuncSeed);
            % Store results
            allResults{ithNet}{ithTest}.replaytrajectory = replaytrajectory{day}{decodeEpoch};
            allResults{ithNet}{ithTest}.ripple.starttime = trialResults.events(:,1).*modelParam.dt;
            allResults{ithNet}{ithTest}.ripple.endtime = trialResults.events(:,2).*modelParam.dt;
        else
            allResults{ithNet}{ithTest}.replaytrajectory = [];
        end

    end
end
end