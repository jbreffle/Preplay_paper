function [PFSpikeTimes, preplaySpikeTimes, rngStruct, avg_mat, allResults, PFresults] = run_param_sim(modelParam, simParam, ithParamSet, mainSeed, varargin)
% run_param_sim
%
% Simulate a set of networks at a given parameter point:
% This function replaced the old function parallelize_parameter_tests_2_PF
%
% Inputs:
%
%
% Outputs:
%
%

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'simParam',	@isstruct)
addRequired(inputObj, 'ithParamSet',	@isnumeric)
addRequired(inputObj, 'mainSeed',	@isnumeric)
parse(inputObj, modelParam, simParam, ithParamSet, mainSeed, varargin{:});
p = inputObj.Results;


%% Initialization

% Set up parameter values for current parameter set
for i = 1:size(simParam.variedParam, 2)
    modelParam.(simParam.variedParam(i).name) = simParam.parameterSets_vec(i,ithParamSet);
end

% Update any parameters that are dependent on a varied parameter
modelParam = set_depedent_parameters(modelParam);

%Run network initialization code
resp_mat = zeros(simParam.nNets, 4);
allResults = cell(1, simParam.nNets) ;
PFresults = cell(1, simParam.nNets) ;

% New, to save all spikes
PFSpikeTimes = cell(simParam.nNets, modelParam.nEnvironments, modelParam.nTrials_PF, modelParam.n);
preplaySpikeTimes = cell(simParam.nNets, simParam.nTrials_preplay, modelParam.n);
rngStruct = struct;
rngStruct.mainSeed = mainSeed;

%% Run sim
for ithNet = 1:simParam.nNets
    % disp(['Starting net ', num2str(ithNet), ' of paramset ', num2str(ithParamSet), ' at ', datestr(datetime(now, 'ConvertFrom', 'datenum'))])

    rngStruct.net(ithNet).netSeed = mainSeed*ithNet;
    network = create_network(modelParam, 'seed', rngStruct.net(ithNet).netSeed);

    rngStruct.net(ithNet).E_indices = network.E_indices;
    rngStruct.net(ithNet).parameterSets_vec = simParam.parameterSets_vec(:,ithParamSet);
    PFresults{ithNet}{1}.E_indices = network.E_indices; % only collect E_indices to verify that create_clusters with netSeed is successful at recreating networks posthoc
    PFresults{ithNet}{1}.netSeed = ithNet;
    PFresults{ithNet}{1}.variedParamVals = simParam.parameterSets_vec(:,ithParamSet);

    mat = zeros(simParam.nTrials_preplay,4);

    %% Run PF simulation

    %opV = zeros(parameters.n, numel(modelParam.t_PF), modelParam.nEnvironments, modelParam.nTrials_PF); % Voltage from all PF sims
    opS = false(modelParam.n, numel(modelParam.t_PF), modelParam.nEnvironments, modelParam.nTrials_PF); % Spikes from all PF sims
    for ithEnv = 1:modelParam.nEnvironments
        for ithTrial = 1:modelParam.nTrials_PF

            % Set up for ithTrial PF simulation
            isPFSim = true;
            rngStruct.net(ithNet).PF_G(ithTrial) = ithTrial*ithNet*mainSeed;
            trialGinSeed = rngStruct.net(ithNet).PF_G(ithTrial);

            % This is needed to make sure V_m0 starts from the same rng state as the simulations in the preprint
            % Needed due to switching to the randnet_calculator_memOpt_wGin() function
            rng(trialGinSeed, 'twister')
            nTimeSteps = ( modelParam.t_max_PF/modelParam.dt)+1;
            randn(modelParam.n, 1);
            for i = 2:nTimeSteps
                rand(modelParam.n, 1); rand(modelParam.n, 1); rand(modelParam.n, 1); rand(modelParam.n, 1);
            end
            V_m0 = modelParam.vmInitMean + modelParam.vmInitSTD*randn([modelParam.n,1]);
            G_sra0 = zeros(modelParam.n, 1);

            % PF Simulation
            spikeMat_PF = randnet_calculator_memOpt_wGin(modelParam, network, ithEnv, isPFSim, seed=trialGinSeed, V_m0=V_m0, G_sra0=G_sra0);

            opS(:,:,ithEnv,ithTrial) = spikeMat_PF;

            % Collect spike times to save
            for ithCell=1:modelParam.n
                spikeInd_dt = find(spikeMat_PF(ithCell,:));
                PFSpikeTimes{ithNet, ithEnv, ithTrial, ithCell} = spikeInd_dt.*modelParam.dt;
            end

        end
    end

    %% Analyze PF simulation
    if simParam.run_analysis_within_sim>0
        linfields = calculate_linfields(opS, modelParam);
        [PFmat, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);
        PFresults{ithNet}{1}.linfields = PFmat;

        allScores = zeros(1, modelParam.nEnvironments);
        if modelParam.PFscoreFlag % Calculating the score is computationally expensive (due to fitting gaussian curves)
            for ithEnv = 1:modelParam.nEnvironments
                allScores(ithEnv) = calculate_linfieldsScore(linfields, modelParam, network, trajectory=ithEnv);
                disp(['Env: ', num2str(i), ', Score: ', num2str(allScores(i))])
            end
            disp(['Mean score: ', num2str(mean(allScores))])
        end
    else
        PFresults{ithNet}{1}.linfields = [];
    end


    clear G_in_PFs opV opS


    %% Preplay simulations

    for ithTest = 1:simParam.nTrials_preplay
        %% Run a preplay simulation
        rngStruct.net(ithNet).preplay_G(ithTest) = ithTest*mainSeed;
        trialGinSeed = rngStruct.net(ithNet).preplay_G(ithTest);
        isPFSim = false;
        preplayContext = 1; % Which environments context cue to use for preplay sim

        %Run model
        V_m0 = repmat(modelParam.V_reset, modelParam.n, 1); % set all neurons to baseline reset membrane potential
        G_sra0 = zeros(modelParam.n, 1);
        spikeMat = randnet_calculator_memOpt_wGin(modelParam, network, preplayContext, isPFSim, seed=trialGinSeed, V_m0=V_m0, G_sra0=G_sra0, initialize=false);

        % TODO: vectorize?
        % Collect spike times to save
        for ithCell=1:modelParam.n
            spikeInd_dt = find( spikeMat(ithCell,:) );
            preplaySpikeTimes{ithNet, ithTest, ithCell} = spikeInd_dt.*modelParam.dt;
        end

        %% Preplay analysis

        if simParam.run_analysis_within_sim>0
            %% Basic preplay statistics
            E_spikes_V_m = spikeMat(network.E_indices,:); clear spikeMat

            % detect events and compute outputs
            [trialResults] = detect_PBE(E_spikes_V_m, modelParam);
            if ithTest == 1 % append trialResults struct to network results struct
                network_spike_sequences = trialResults;
            else
                network_spike_sequences = [network_spike_sequences, trialResults];
            end

            % Calculate mean n spikes within events
            overallmeanspikes = 0;
            for ithEvent = 1:size(trialResults.events, 1)
                eventSpikes = E_spikes_V_m(:,trialResults.events(ithEvent,1):trialResults.events(ithEvent,2));
                nSpikes = sum(eventSpikes, 2);
                meannSpikes = sum(nSpikes)/sum(nSpikes>0);
                overallmeanspikes = overallmeanspikes+ meannSpikes;
            end
            overallmeanspikes = overallmeanspikes./size(trialResults.events, 1);

            % Overall simulation statistics
            allResults{ithNet}{ithTest}.ithInit = ithTest;
            allResults{ithNet}{ithTest}.numEvents = numel(network_spike_sequences(ithTest).event_lengths); % number of detected events
            allResults{ithNet}{ithTest}.fracFire =  mean(sum(E_spikes_V_m, 2)>0); % Fraction of cells that fire at all during simulation
            allResults{ithNet}{ithTest}.frac_participation = mean([network_spike_sequences(ithTest).frac_spike{:}]); % mean fraction of cells firing per event
            allResults{ithNet}{ithTest}.meanRate = mean(sum(E_spikes_V_m, 2)/modelParam.t_max_preplay); % mean over cells' average firing rate
            allResults{ithNet}{ithTest}.stdRate = std(sum(E_spikes_V_m, 2)/modelParam.t_max_preplay); % STD over cells' average firing rate
            allResults{ithNet}{ithTest}.meanCellnEventSpikes = overallmeanspikes; % mean n event spikes over cells that spike

            % Stats for each detected event
            allResults{ithNet}{ithTest}.eventLength = network_spike_sequences(ithTest).event_lengths; % duration in seconds of all detected events
            allResults{ithNet}{ithTest}.eventParticipation = [network_spike_sequences(ithTest).frac_spike{:}]; % fraction of cells that fired in each event

            % Save ranks_vec
            allResults{ithNet}{ithTest}.ranksVec = network_spike_sequences(ithTest).ranks_vec;


            %% Calculate and save preplay Bayesian decoding
            if ~isempty(trialResults.events) % only decode if there are events
                % 1) save SJlab formated files to tmp dir
                [filepath,name] = fileparts( tempname([pwd, '/results/tmp/']) );
                animalprefix = name; save_path = [filepath, filesep, animalprefix, '_direct'];
                mkdir(save_path)

                day = 1; eprun = 2; tet = 1;
                tetinfo{day}{eprun}{tet}.numcells = modelParam.n_E;
                if numel(linfields{day})<eprun
                    linfields{day}{eprun} = linfields{day}{1};
                end
                if numel(linfields{day}{eprun}{:})==modelParam.n
                    linfields{day}{eprun}{tet}(network.I_indices) = [];
                end

                ripple{day}{eprun}.starttime = trialResults.events(:,1).*modelParam.dt;
                ripple{day}{eprun}.endtime = trialResults.events(:,2).*modelParam.dt;
                for cellid = 1:numel(network.E_indices)
                    % spikes{day}{eprun}{tet}{cellid}.data(:,1) = V_m(network.E_indices(cellid),:)>= parameters.V_th; % spike times
                    spikes{day}{eprun}{tet}{cellid}.data(:,1) = find(E_spikes_V_m(cellid,:)).*modelParam.dt; % spike times
                end
                clear E_spikes_V_m

                save([filepath, filesep, animalprefix, '_direct/', animalprefix, 'tetinfo.mat'], 'tetinfo')
                save([filepath, filesep, animalprefix, '_direct/', animalprefix, 'linfields01.mat'], 'linfields')
                save([filepath, filesep, animalprefix, '_direct/', animalprefix, 'rippletime01.mat'], 'ripple')
                save([filepath, filesep, animalprefix, '_direct/', animalprefix, 'spikes01.mat'], 'spikes')
                clear tetinfo linfields ripple spikes

                % 2) Save results:
                day=1; ep=2; savedata = true; figopt = 0;
                rngStruct.net(ithNet).decodeFunc(ithTest) = ithTest*ithNet*1000*mainSeed;
                decodeFuncSeed = rngStruct.net(ithNet).decodeFunc(ithTest);
                if simParam.run_analysis_within_sim>1
                    warning off
                    decodingResults = decode_events(modelParam, animalprefix,[filepath, filesep], savedata, day=day, epoch=ep, figopt=figopt, decodeFuncSeed=decodeFuncSeed);
                    warning on
                    allResults{ithNet}{ithTest}.replaytrajectory = decodingResults{day}{ep};
                end
                % 3) Delete tmp dir files:
                rmdir(save_path, 's')
            else
                allResults{ithNet}{ithTest}.replaytrajectory = [];
            end

            % Basic output statistics
            mat(ithTest,1) = allResults{ithNet}{ithTest}.frac_participation; % fraction of spiking neurons over identified events
            mat(ithTest,2) = allResults{ithNet}{ithTest}.meanRate ; %average firing rate
            mat(ithTest,3) = mean(allResults{ithNet}{ithTest}.eventLength); % Average event length
            mat(ithTest,4) = allResults{ithNet}{ithTest}.numEvents; % Number of events
        else
            mat(ithTest,:) = [];
            allResults{ithNet}{ithTest} = [];
        end


    end % trial loop
    mat(isnan(mat)) = 0;
    resp_mat(ithNet,:) = sum(mat,1) ./ sum(mat > 0,1); %Only averaging those that did successfully produce data

end % Network loop
resp_mat(isnan(resp_mat)) = 0;
avg_mat = sum(resp_mat,1)./sum(resp_mat > 0,1); %Only averaging those that had results
% TODO: delete avg_mat, resp_mat, mat

%disp(datetime(now, 'ConvertFrom', 'datenum'))
%disp(['Parameter set ', num2str(ithParamSet), '/', num2str(size(parameterSets_vec, 2)), ' complete'])

end % Function