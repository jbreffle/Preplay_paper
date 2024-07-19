%% sim_randnet.m
%
% Simulates both PFs and preplay for a single network at a single
% parameter point.

plotFig2 = true; % Set to true to generate Figure 2a-b of the manuscript


%% Automatically get data location

addpath(fullfile("src"))
Config = utils.setConfig;

[s, simParam.sim_git_hash_string] = system('git rev-parse HEAD');

myPlotSettings

% Set up name and dir for saving results (if simParam.saveFlag=true)
simFileName = ['randnet_', datestr(now,'yyyy-mm-ddTHH-MM')];
simFilePath = fullfile(Config.simDataPath, 'sim-randnet', simFileName);


%% Simulation options

simParam.saveFlag = false;      % If true, save simulation results
simParam.saveGvar = false;      % If true, save all conductances from the preplay
simParam.saveVm = false;        % If true, save all Vm from PF and preplay sims

simParam.plotResults = true;    % If true, plot simulation results
simParam.calcPFScore = false;   % If true, calculate PF score

simParam.runPFSim = true;       % If true, run PF simulation
simParam.runPreplaySim = true;	% If true, run preplay simulation
simParam.runDecode = true;     % If true, run decoding using the PF and preplay simulations

% Simulation options that are passed to the modelParam struct
simParam.nNets          = 1;    % Number of networks to simulate at each parameter point
simParam.nTrials_preplay= 1;    % Number of preplay trials to simulate per network
simParam.t_max_preplay  = 10;   % Trial duration (s) for each preplay trial
simParam.nTrials_PF     = 5;    % Number of PF trials to simulate (each run across the track)
simParam.t_max_PF       = 2;    % Duration (s) for each PF trial
simParam.envIDs         = [1];	% Environmnent IDs to simulate for PFs (default is [1:4] if this is empty or non-existent)

% New
% simParam.nParPool = 4;
simParam.dispFlag_decode = 0;

MSPlots = true; % If true, plot panels for manuscript
%figSettings = 'manuscript';
%assert(any(strcmp(figSettings, {'standard', 'manuscript'})))
netSeed = 1;

% Override settings to generate the manuscript figures
if plotFig2
    disp('Plotting manuscript figures')
    netSeed = 1;
    simParam.nTrials_PF = 1;
    simParam.envIDs = 1:14;
    simParam.runPreplaySim = false;
    simParam.plotResults = false;
    MSPlots = true;
end


%% Initialize parameters

simParam.paramSet = '7_22-v2' ;	% 'default', '7_22', '10_18', 'ISregMin', 'PFtuned1', 'PFtuned2', 'PFtuned3'
disp(['Using ', simParam.paramSet, ' parameters'])
modelParam = initialize_parameters(simParam.paramSet, simParam);


%% Create network

network = create_network(modelParam, 'seed', netSeed);


%% Place field simulation:
tic
if simParam.runPFSim

    opV_PF = zeros(modelParam.n, numel(modelParam.t_PF), modelParam.nEnvironments, modelParam.nTrials_PF);	% Voltage from all PF sims
    opGin_PF = zeros(modelParam.n, numel(modelParam.t_PF), modelParam.nEnvironments, modelParam.nTrials_PF); % Input conductances for all PF sims
    opS_PF = false(modelParam.n, numel(modelParam.t_PF), modelParam.nEnvironments, modelParam.nTrials_PF);	% Voltage from all PF sims
    for ithEnv = 1:modelParam.nEnvironments
        for ithTrial = 1:modelParam.nTrials_PF
            % Set up for ithTrial PF simulation
            isPFSim = true;
            trialGinSeed = ithEnv*100+ithTrial;
            G_in_PFs = create_G_in(modelParam, network, ithEnv, isPFSim, seed=trialGinSeed, runTests=true);
            V_m0 = modelParam.vmInitMean + modelParam.vmInitSTD*randn([modelParam.n,1]);
            G_sra0 = zeros(modelParam.n, 1);

            % PF Simulation
            [spikeMat_PF, V_m_PF, ~, ~, ~, ~, ~] = randnet_calculator(modelParam, network, modelParam.t_max_PF, G_in=G_in_PFs, V_m0=V_m0, G_sra0=G_sra0);
            opV_PF(:,:,ithEnv,ithTrial) = V_m_PF;
            opS_PF(:,:,ithEnv,ithTrial) = spikeMat_PF;
            opGin_PF(:,:,ithEnv,ithTrial) = G_in_PFs;
        end
    end

    % Calculate place fields, following Jadhav lab structure
    linfields = calculate_linfields(opS_PF, modelParam);

    % Calculate PF score
    if simParam.calcPFScore
        allScores = zeros(1, modelParam.nEnvironments);
        for ithEnv = 1:modelParam.nEnvironments
            allScores(ithEnv) = calculate_linfieldsScore(linfields, modelParam, network, trajectory=ithEnv);
            disp(['Env: ', num2str(ithEnv), ', Score: ', num2str(allScores(ithEnv))])
        end
        disp(['Mean score: ', num2str(mean(allScores))])
    end

    % Plot PF results
    if simParam.plotResults
        % Plot example data
        plots.randnetPlaceFields(modelParam, network, opS_PF, opV_PF, opGin_PF, linfields, ithEnv=1, ithTrial=1)
        % Plot the place fields
        plots.linfieldsPlaceFields(linfields, network, modelParam)
    end
    if MSPlots
        plots.randnetPlaceFieldsManuscript(modelParam, network, opS_PF, opV_PF, opGin_PF, linfields, ithEnv=1, ithTrial=1)
    end

    plotMSRawPF = [modelParam.nEnvironments>=14];
    % TODO: clean up and move into function
    if plotMSRawPF
        assert(modelParam.nEnvironments>4)
        assert(simParam.nTrials_PF==1)
        myPlotSettings(width=2.25, height=1.5)
        MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
        [~, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);
        for ithPlotEnv = 14%1:modelParam.nEnvironments
            for ithPlotTrial = 1:modelParam.nTrials_PF
                figure;
                Espikes = logical(opS_PF(network.E_indices(PFpeaksSequence{ithPlotEnv}),:,ithPlotEnv,ithPlotTrial));
                %rpermIcells = 1:numel(network.I_indices); %randperm(numel(network.I_indices)); % Randomly permute cells, to prove they are not sorted
                %plotSpikeRaster( logical( [ opS_PF(network.E_indices(PFpeaksSequence{ithPlotEnv}),:,ithPlotEnv,ithPlotTrial); opS_PF(network.I_indices(rpermIcells),:,ithPlotEnv,ithPlotTrial) ]), 'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
                plotSpikeRaster(Espikes, 'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat);
                xlabel('Time (s)'); ylabel('Cell (sorted)');
                title(['Env ', num2str(ithPlotEnv), ', Trial ', num2str(ithPlotTrial)])
                ylim([1-2, find( any(Espikes, 2) , 1 , 'last' )+2])

                myPlotSettings(width=1.5, height=1.5)
                sortedCellIndsToPlot = [91, 100, 109, 117, 124, 140]; % Choose cells here
                cellInds = network.E_indices(PFpeaksSequence{ithPlotEnv}(sortedCellIndsToPlot));
                for cellInd_i = cellInds
                    cellVm = opV_PF(cellInd_i,:,ithPlotEnv,ithPlotTrial);
                    cellVm(cellVm==modelParam.V_reset) = 40e-3;
                    figure; plot(modelParam.t_PF, 1000*cellVm)
                    yAxisRanges = [-59 -45]; ylim(yAxisRanges)
                    box off; ylabel('V_m (mV)'); xlabel('Time (s)');
                    title(['Env ', num2str(ithPlotEnv), ', Trial ', num2str(ithPlotTrial), ', cell ', num2str(cellInd_i)])
                end
            end
        end
        myPlotSettings
    end

    PFRuntime = toc;
    disp(['PF sim runtime: ', num2str(PFRuntime), ' s'])
else
    disp('Skipping PF simulation')
end


%% Preplay simulation:
tic
if simParam.runPreplaySim

    opV_preplay = zeros(modelParam.n, modelParam.t_steps_preplay, modelParam.nTrials_preplay); % Voltage from all preplay sims
    opS_preplay = false(modelParam.n, modelParam.t_steps_preplay, modelParam.nTrials_preplay); % Spikes from all preplay sims
    G_var = struct;
    for ithTest = 1:modelParam.nTrials_preplay
        % Set up for ithTest preplay simulation
        preplayContextEnv = 1; % Which environments context cue to use for preplay sim, unless modelParam.useSleepContext==true
        isPFSim = false;
        trialGinSeed = ithTest;
        G_var(ithTest).G_in = create_G_in(modelParam, network, preplayContextEnv, isPFSim, seed=trialGinSeed, runTests=true);
        V_m0 = modelParam.vmInitMean + modelParam.vmInitSTD*randn([modelParam.n,1]);
        G_sra0 = zeros(modelParam.n, 1);

        % Run model
        [spikeMat, V_m, G_var(ithTest).G_sra, G_var(ithTest).G_syn_E_E, G_var(ithTest).G_syn_I_E, G_var(ithTest).G_syn_E_I, G_var(ithTest).G_syn_I_I] = ...
            randnet_calculator(modelParam, network, modelParam.t_max_preplay, G_in=G_var(ithTest).G_in, V_m0=V_m0, G_sra0=G_sra0);
        opV_preplay(:,:,ithTest) = V_m;
        opS_preplay(:,:,ithTest) = spikeMat;

        % Detect PBEs
        trialResults = detect_PBE(spikeMat(network.E_indices,:), modelParam, plotAll=false);
        if ithTest == 1 %  Append trialResults struct to network results struct
            network_spike_sequences = trialResults;
        else
            network_spike_sequences = [network_spike_sequences, trialResults];
        end

        % Plot preplay results
        if simParam.plotResults
            plots.randnetPreplay(modelParam, network, V_m, G_var(ithTest).G_in, network_spike_sequences(ithTest), sortMethod=2)
        end

    end % Trial loop
    preplayRuntime = toc;
    disp(['Preplay sim runtime: ', num2str(preplayRuntime), ' s'])
else
    disp('Skipping preplay simulation')
end


%% Run decode
if simParam.runDecode
    if exist('trialResults') && ~isempty(trialResults.events) && any(trialResults.event_lengths>=(modelParam.minEventDur*1e-3))
        modelParam.shuffleIterations = 100;
        tic
        replaytrajectory = struct;
        for ithTestToDecode = 1%:simParam.nTrials_preplay
            trialDecodeStruct = run_net_decode(modelParam, network, trialResults, linfields, squeeze(opS_preplay(:,:,ithTestToDecode)), figOpt=2);
            replaytrajectory = trialDecodeStruct;
            plots.randnetDecode(trialDecodeStruct)
        end
        decodeRuntime = toc;
        disp(['Decode runtime: ', num2str(decodeRuntime), ' s'])
    else
        disp('No events to decode')
    end
end


%%

if simParam.runDecode && exist('trialResults') && ~isempty(trialResults.events) && any(trialResults.event_lengths>=(modelParam.minEventDur*1e-3))
    [PFmat, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);
    rippLength = (trialResults.events(:,2)- trialResults.events(:,1))*modelParam.dt;
    goodRippInd =  find(rippLength> (modelParam.minEventDur/1000));
    ithEvent = 1;
    ithRipp = goodRippInd(ithEvent);
    PFmatAllTraj = [PFmat{:}];
    PFmat_E = PFmatAllTraj(network.E_indices,:);
    PFpeaks = max(PFmat_E, [], 2);
    EcellsSpiked = opS_preplay(network.E_indices,trialResults.events(ithRipp,1):trialResults.events(ithRipp,2));
    ECellValid = PFpeaks>=modelParam.minPeakRate;
    rasterCells = find(any(EcellsSpiked & ECellValid, 2))';
    decodeCells1 = replaytrajectory{1}{2}.activecell{ithEvent}';
    decodeCells2 = replaytrajectory{1}{2}.activecellidx{ithEvent}(:,2)';
    assert( isequal(rasterCells, decodeCells2))
    disp('decode cells validated')
end


%% Save results

if simParam.saveFlag
    if ~isfolder(simFilePath)
        mkdir(simFilePath);
    end

    % Save network structure
    save(strcat(simFilePath,'/network.mat'),'network');

    % Save parameter files
    save(strcat(simFilePath,'/allParams.mat'), 'parameters', 'pfsim', 'simParam');

    % Save simulation data
    save(strcat(simFilePath,'/spike_mats.mat'), 'opS_preplay', 'opS_PF') %,'-v7.3')
    save(strcat(simFilePath,'/network_spike_sequences.mat'),'network_spike_sequences') %,'-v7.3')

    % Save decode results if they were produced
    if exist('replaytrajectory', 1)
        save(strcat(simFilePath,'/decodeResult.mat'), 'replaytrajectory');
    end

    % Save Vm and conductances if requested
    if simParam.saveVm
        save(strcat(simFilePath,'/Vm_mats.mat'), 'opV_preplay', 'opV_PF') %,'-v7.3')
    end
    if simParam.saveGvar
        save(strcat(simFilePath,'/G_var.mat'),'G_var', 'opGin_PF') %,'-v7.3')
    end

    disp(['Results saved at: ', newline, simFilePath])
else
    disp('Results not saved')
end

