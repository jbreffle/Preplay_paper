% testRandnetCalculator.m
%
% Adapted from simulateNetwork.m. Tests that all randnet_calculators produce
% identical spike results.
%


%% Set up

addpath("../src")
Config = utils.setConfig();

[s,simParam.sim_git_hash_string] = system('git rev-parse HEAD');

myPlotSettings


%% Simulation options

simParam.saveFlag = false;      % If true, save simulation results
simParam.saveGvar = false;      % If true, save all conductances from the preplay
simParam.saveVm = false;        % If true, save all Vm from PF and preplay sims

simParam.plotResults = true;    % If true, plot simulation results

simParam.runPFSim = true;       % If true, run PF simulation
simParam.runPreplaySim = true;	% If true, run preplay simulation
simParam.runDecode = true;     % If true, run decoding using the PF and preplay simulations

% Simulation options that are passed to the modelParam struct
simParam.calcPFScore = false;   % If true, calculate PF score
simParam.nNets          = 1;    % Number of networks to simulate at each parameter point
simParam.nTrials_preplay= 1;    % Number of preplay trials to simulate per network
simParam.t_max_preplay  = 2;   % Trial duration (s) for each preplay trial
simParam.nTrials_PF     = 5;    % Number of PF trials to simulate (each run across the track)
simParam.t_max_PF       = 2;    % Duration (s) for each PF trial

% New
% simParam.nParPool = 4;
simParam.dispFlag_decode = 0;


%% Initialize parameters

paramSet = '7_22' ;	% 'default', '7_22', '10_18', 'ISregMin', 'PFtuned'
disp(['Using ', paramSet, ' parameters'])
modelParam = initialize_parameters(paramSet, simParam);


%% Create network
netSeed = 1;
network = create_network(modelParam, 'seed', netSeed);


%% Preplay simulation:
tic

for ithTest = 1:modelParam.nTrials_preplay

    % Set up for ithTest preplay simulation
    preplayContextEnv = 1; % Which environments context cue to use for preplay sim
    isPFSim = false;
    trialGinSeed = ithTest;
    G_in = create_G_in(modelParam, network, preplayContextEnv, isPFSim, seed=trialGinSeed, runTests=true);
    V_m0 = modelParam.vmInitMean + modelParam.vmInitSTD*randn([modelParam.n,1]);
    G_sra0 = zeros(modelParam.n, 1);

    %% Test old vs updated
    tic;
    [spikeMat_main, V_m_main] = randnet_calculator(modelParam, network, modelParam.t_max_preplay, G_in=G_in, V_m0=V_m0, G_sra0=G_sra0);
    toc

    V_m = zeros(modelParam.n, modelParam.t_steps_preplay); V_m(:,1) = V_m0;
    modelParam.G_in=G_in; trialSeed = randi(10^6);
    tic;
    [V_m_temp] = randnet_calculator_temp(modelParam, trialSeed, network, V_m, modelParam.t_max_preplay);
    toc

    spikes_main = logical( V_m_main>modelParam.V_th);
    spikes_temp = logical( V_m_temp>modelParam.V_th);
    assert(isequal(V_m_main, V_m_temp))
    assert(isequal(spikes_main, spikes_temp))

    %% Test updated vs memOpt
    tic;
    [spikes_memOpt] = randnet_calculator_memOpt(modelParam, network, modelParam.t_max_preplay, G_in=G_in, V_m0=V_m0, G_sra0=G_sra0);
    toc
    assert(isequal(spikes_main, spikes_memOpt))
    assert(isequal(spikes_main, spikeMat_main))

    % Speed comparison
    tic
    for ii = 1:5; [spikeMat_main] = randnet_calculator(modelParam, network, modelParam.t_max_preplay, G_in=G_in, V_m0=V_m0, G_sra0=G_sra0); end
    disp(['Standard 5x runtime: ', num2str(toc)])
    tic
    for ii = 1:5; [spikes_memOpt] = randnet_calculator_memOpt(modelParam, network, modelParam.t_max_preplay, G_in=G_in, V_m0=V_m0, G_sra0=G_sra0); end
    disp(['Mem-opt  5x runtime: ', num2str(toc)])

    %% Test memOpt vs memOpt_wGin
    tic
    [spikes_memOpt] = randnet_calculator_memOpt(modelParam, network, modelParam.t_max_preplay, G_in=G_in, V_m0=V_m0, G_sra0=G_sra0);
    toc
    tic
    [spikes_memOpt_wGin] = randnet_calculator_memOpt_wGin(modelParam, network, preplayContextEnv, isPFSim, seed=trialGinSeed, V_m0=V_m0, G_sra0=G_sra0);
    toc
    assert(isequal(spikes_memOpt_wGin, spikes_memOpt))

end % Trial loop
disp('Tests successful!')
