function modelParam = initialize_parameters(paramSet, simParam, varargin)
% initialize_parameters.m
%
% Creates the baseline parameter structure, then original parameters are
% optionally overridden based on paramSet.
%
% Input:
%   paramSet: String choice of parameter set {'default', '7_22', or '10_18'}
%   simParam: simulation parameter structure (this dependency might be
%   remove later, by removing the corresponding fields from modelParam)
%
% Output:
%   The parameters structure
%

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'paramSet',	@ischar)
addRequired(inputObj, 'simParam',	@isstruct)
parse(inputObj, paramSet, simParam, varargin{:});
p = inputObj.Results;

%% Original model parameters

% Simulation timestep
modelParam.dt = 0.1*10^(-3); %timestep (s)

% Network structure parameters
modelParam.n = 500;         % Number of neurons
modelParam.clusters = 10;	% Number of clusters
modelParam.mnc = 1.5;         % Mean number of clusters each neuron is a member of

% Basic model parameters
modelParam.tau_syn_E = 10*10^(-3);  % Exc. synaptic decay time constant (s), ~10 ms from DOI:10.1126/science.aaf1836
modelParam.tau_syn_I = 3*10^(-3);	% Inh. synaptic decay time constant (s), ~1.2-8 ms from DOI:10.1073/pnas.192233099
modelParam.tau_stdp = 5*10^(-3);	% STDP time constant (s)

modelParam.E_K = -80*10^(-3);       % Potassium reversal potential (V)
modelParam.E_L = -70*10^(-3);       % Leak reversal potential (V), -60 to -70 mV range
modelParam.G_L = 25*10^(-9);        % Leak conductance (S), 10 to 30 nS range
modelParam.C_m = 0.4*10^(-9);       % Total membrane capacitance (F) %Huge range from 0.1 - 100 pF
modelParam.V_th = -50*10^(-3);      % Spike threshold membrane potential (V)
modelParam.V_reset = -70*10^(-3);   % Spike reset membrane potential (V)
modelParam.V_syn_E = 0;             % Synaptic reversal potential (excitatory)
modelParam.V_syn_I = -70*10^(-3);	% Synaptic reversal potential (inhibitory)

% Recurrent connection strengths
modelParam.del_G_syn_E_E = 750*10^(-12);	% E to E synaptic conductance step following spike (S)
modelParam.del_G_syn_I_I = 0;               % I to I synaptic conductance step following spike (S)
modelParam.del_G_syn_E_I = 500*10^(-12);    % synaptic conductance step following spike (S)
modelParam.del_G_syn_I_E = nan;             % synaptic conductance step following spike (S)

% Spike rate adaptation (SRA) parameters
modelParam.del_G_sra = 330e-09;     % SRA conductance step following spike %ranges from 1-200 *10^(-9) (S)
modelParam.tau_sra = 30*10^(-3);	% SRA time constant (s)

% Poisson input parameters:
modelParam.rG = 1000;               % Poisson input spiking rate
modelParam.Win_mean = 73 *10^-12;   % Poisson input weight mean
modelParam.Win_var = (5e-12)^2;     % Poisson input weight variance

% Network connection parameters
modelParam.conn_prob = 0.08;    % Global E-to-E connection probability
modelParam.p_E = 0.75;          % Fraction of cells that are excitatory
modelParam.include_all = 2;     % If>1, neurons not in a cluster are added to one (different methods)
modelParam.global_inhib = true;    % If 1, I-cells are not clustered and have connection probability p_I
modelParam.p_I = 0.25;          % I-to-I, E-to-I, and I-to-E connection probability

% New parameters
modelParam.clusterCorrs = 1;        % If 1, cluster-based correlations are added  to the input weights
modelParam.inputBiasSigma = 25;     % 1/inputBiasSigma is the biased added if clusterCorrs=1
modelParam.PFcontextFrac = 0.1;     % Fraction of input that is from the context signal for PF simulations

% Initial mean and STD of membrane voltages
modelParam.vmInitMean  = -52.5e-3;
modelParam.vmInitSTD   = 1e-3;


%% PF specific parameters
if isfield(simParam, 'envIDs') && ~isempty(simParam.envIDs)
    modelParam.envIDs = simParam.envIDs;
else
    modelParam.envIDs = [1, 2, 3, 4] ; % odd IDs are rightward traversals, even IDs are rightward traversals of the ID-1 track
end
modelParam.trackWidth = 1;           % Maximum value (m)
modelParam.spatialBin = 2/100;       % 2 cm bins for localizing position
modelParam.linFieldGaussSD = 0.04;   % Standard deviation of PF gaussian kernel
modelParam.winScale = 5;             % Window is 5x the STDev


%% Analysis parameters (split into analysisParam structure?)

% PBE detection parameters
modelParam.PBE_min_Hz = 0.5; % Minimum population mean rate during PBE, to exclude pathologically silent parameter sets
modelParam.PBE_zscore = 1.0; % minimum stds above mean rate to detect PBE
modelParam.PBE_min_dur = 30 * (1/1000); % Minimum duration of a PBE (modelParam.minEventDur is minimum for decode)
modelParam.PBE_window =  15 * (1/1000) *(1/modelParam.dt); % Width of gaussian kernel used to calculate mean pop activity rate
modelParam.PBE_max_combine = 10 * (1/1000);     % Combine adjacent PBEs separaeted by less than this duration
modelParam.PBE_useMeanCrossing = false;         % If true, use the mean-crossing points before/after event as the start/end

% PF Gaussian fit parameters
modelParam.gaussFOLower = [10, 0, 0]; % [peak amplitude, position of peak on track, standard deviation of peak]
modelParam.gaussFOUpper = [30, modelParam.trackWidth*100, sqrt(modelParam.trackWidth*100)];
modelParam.peakTarget = 15; % Target peak rate for linfieldsScore (Hz)

% Decoding parameters
modelParam.cellcountthresh = 5;     % Minimum number of participating cells for a candidate event
modelParam.shuffleIterations = 500;	% Number of time-bin shuffles for each event
modelParam.useLogSumDecode = true;	% Use sum of log probabilities for decoding events
modelParam.tBinSz = 10;            % ms, time bin size for decoding
modelParam.minEventDur = 50;       % ms, exclude events shorter than this
modelParam.wellcutoff = 0;         % cm, remove reward-well regions (15cm around start and end); or 0 cm without exclusion
modelParam.minPeakRate = 3;        % Hz, minimum peak rate to include cell as Place Cell
modelParam.normByTraj_decode = true; % Should be true (false is equivalent to old method, which was wrong)
modelParam.downSampleCells = false; % true/false flag to use downsampling
modelParam.downSampledFrac = 0.1;  % fraction of cells to use in decode, if downSampleCells==true


%% simParam dependent settings
modelParam.nNets = simParam.nNets;                  % Number of networks to simulate at each parameter point
modelParam.nTrials_preplay = simParam.nTrials_preplay;	% Number of preplay trials to simulate per network
modelParam.t_max_preplay = simParam.t_max_preplay;	% Trial duration (s) for each preplay trial
modelParam.PFscoreFlag = simParam.calcPFScore;      % If true, calculate PF score
modelParam.nTrials_PF = simParam.nTrials_PF;        % Number of PF trials to simulate (each run across the track)
modelParam.t_max_PF = simParam.t_max_PF;            % Duration (s) for each PF trial
modelParam.t_PF = 0:modelParam.dt:modelParam.t_max_PF;
modelParam.xPos = 0:modelParam.dt/modelParam.t_max_PF:1;
modelParam.dispFlag_decode = simParam.dispFlag_decode;


%% Standardized parameter set changes
switch paramSet
    case 'default'
        % do nothing
    case '7_22'
        modelParam.del_G_syn_E_E =  220e-12;
        modelParam.del_G_syn_E_I =  400e-12;
        modelParam.Win_mean = 72e-12;

        modelParam.IcueScale_PF = 1;
        modelParam.IcueScale = 0.7500;
        modelParam.inhSRA = 0;

        modelParam.G_L = 10*10^(-9);
        modelParam.del_G_sra = 3.0e-012;
        modelParam.rG = 5000;
    case '7_22-v2'
        % Identical to 7_22, but with useSleepContext = true;
        modelParam.del_G_syn_E_E =  220e-12;
        modelParam.del_G_syn_E_I =  400e-12;
        modelParam.Win_mean = 72e-12;

        modelParam.IcueScale_PF = 1;
        modelParam.IcueScale = 0.7500;
        modelParam.inhSRA = 0;

        modelParam.G_L = 10*10^(-9);
        modelParam.del_G_sra = 3.0e-012;
        modelParam.rG = 5000;

        modelParam.useSleepContext = true; % Only change to 7_22

    case '10_18'
        modelParam.del_G_syn_E_E = 150*10^(-12);
        modelParam.del_G_syn_E_I = 90*10^(-12);
        modelParam.Win_mean = 73 *10^-12;

        modelParam.IcueScale_PF = 1.4;
        modelParam.IcueScale = 1.01;
        modelParam.inhSRA = 1;

        modelParam.include_all = 3;
        modelParam.clusters = 8;

        modelParam.G_L = 10*10^(-9);
        modelParam.del_G_sra = 3.0e-012;
        modelParam.rG = 5000;

    case 'PFtuned1'
        % These are the best PF-tuned parameters from early march:
        % They match 'grid_2023-03-06T15-25decode_2023-03-06T16-29'.
        % And 'grid_2023-03-08T09-47decode_2023-03-08T11-34' used the same base arameters
        modelParam.del_G_syn_E_E = 150*10^-12;
        modelParam.del_G_syn_I_E = 225*10^-12;
        modelParam.del_G_syn_E_I = 225*10^-12;
        modelParam.Win_mean = 162 *10^-12 ; % ~double Win_mean
        %modelParam.Win_var = (5e-12)^2;     % Poisson input weight variance

        modelParam.tau_syn_I = 5e-3;
        modelParam.tau_syn_E = 5e-3;
        modelParam.C_m = 0.4*10^(-9) /2;    % Halve C_m

        % These preplay params need tuned:
        modelParam.IcueScale_PF = 1.4;
        modelParam.IcueScale =    1.40;
        modelParam.inhSRA = 1;

        modelParam.include_all = 3;
        modelParam.clusters = 8;
        % modelParam.mnc = 1.5;         % Mean number of clusters each neuron is a member of

        modelParam.G_L = 10*10^(-9);
        modelParam.del_G_sra = 3.0e-012;
        modelParam.rG = 5000;

    case 'PFtuned2'
        % Same as PFtuned1, but with modifications to try to get good preplay
        modelParam.del_G_syn_E_E = 150*10^-12;
        modelParam.del_G_syn_I_E = 225*10^-12;
        modelParam.del_G_syn_E_I = 225*10^-12;
        modelParam.Win_mean = 162 *10^-12 ; % ~double Win_mean
        modelParam.tau_syn_I = 5e-3;
        modelParam.tau_syn_E = 5e-3;
        modelParam.C_m = 0.4*10^(-9) /2;    % Halve C_m
        modelParam.IcueScale_PF = 1.4;
        modelParam.IcueScale =    1.4;
        modelParam.inhSRA = 1;
        modelParam.include_all = 3;
        modelParam.clusters = 8;
        modelParam.G_L = 10*10^(-9);
        modelParam.del_G_sra = 3.0e-012;
        modelParam.rG = 5000;

        % Overwrites of PFtuned1
        modelParam.del_G_syn_E_E = 300*10^-12;
        modelParam.IcueScale_PF = 1.225;
        modelParam.IcueScale    = 1.13;
        modelParam.inhSRA       = 0;
        modelParam.Win_mean     = 150 *10^-12 ; % ~double Win_mean

        % Overwrites of the overwrites
        modelParam.del_G_sra    = 100 * 3.0e-012;


        modelParam.clusters = 15;
        modelParam.mnc = 1.25;

        % Overwrites of the overwrites
        %{
        modelParam.del_G_syn_E_E = 250*10^-12;
        modelParam.IcueScale    = 1.30;
        modelParam.Win_mean     = 162 *10^-12 ; % ~double Win_mean
        %}

        % Overwrites of the overwrites
        %{
        modelParam.del_G_syn_E_E = 250*10^-12;
        modelParam.IcueScale    = 1.30;
        modelParam.Win_mean     = 162 *10^-12 ; % ~double Win_mean
        modelParam.del_G_sra    = 50* 3.0e-012;
        modelParam.inhSRA       = 0;
        modelParam.tau_sra = 100*10^(-3);	% SRA time constant (s)
        %}
        modelParam.vmInitMean  = -55.0e-3;
        disp('initialize_parameters(): PFtuned2 param set under revision')

    case 'PFtuned3'
        % Same as PFtuned1, but with modifications to try to get good preplay
        modelParam.del_G_syn_E_E = 150*10^-12;
        modelParam.del_G_syn_I_E = 225*10^-12;
        modelParam.del_G_syn_E_I = 225*10^-12;
        modelParam.Win_mean = 162 *10^-12 ; % ~double Win_mean
        modelParam.tau_syn_I = 5e-3;
        modelParam.tau_syn_E = 5e-3;
        modelParam.C_m = 0.4*10^(-9) /2;    % Halve C_m
        modelParam.IcueScale_PF = 1.4;
        modelParam.IcueScale =    1.4;
        modelParam.inhSRA = 1;
        modelParam.include_all = 3;
        modelParam.clusters = 8;
        modelParam.G_L = 10*10^(-9);
        modelParam.del_G_sra = 3.0e-012;
        modelParam.rG = 5000;

        % Overwrites of PFtuned1
        modelParam.del_G_syn_E_E = 300*10^-12;
        modelParam.IcueScale_PF = 1.225;
        modelParam.IcueScale    = 1.13;
        modelParam.inhSRA       = 0;
        modelParam.Win_mean     = 150 *10^-12 ; % ~double Win_mean
        modelParam.del_G_sra    = 1 * 3.0e-012;

        modelParam.clusters = 15;
        modelParam.mnc = 1.25;

        % OVerwrites of PFtuned2
        modelParam.del_G_sra    = 1 * 3.0e-012;
        modelParam.del_G_syn_E_E = 250*10^-12;
        % This helps a lot

        % Overwrites of above overwrites
        modelParam.del_G_syn_E_E = 250*10^-12;
        %modelParam.IcueScale    = 1.14;
        %modelParam.clusters = 10;
        %modelParam.mnc = 1.5;
        %modelParam.del_G_sra    = 10 * 3.0e-012;

        modelParam.del_G_syn_I_E = 200*10^-12;
        modelParam.del_G_syn_E_I = nan;
        modelParam.IcueScale_PF = 1.30;
        modelParam.IcueScale    = 1.2;

        modelParam.mnc = 1.75;
        modelParam.del_G_syn_E_E = 290*10^-12;

        % Overwrites of the overwrites
        %{
        modelParam.del_G_syn_E_E = 250*10^-12;
        modelParam.IcueScale    = 1.30;
        modelParam.Win_mean     = 162 *10^-12 ; % ~double Win_mean
        %}

        % Overwrites of the overwrites
        %{
        modelParam.del_G_syn_E_E = 250*10^-12;
        modelParam.IcueScale    = 1.30;
        modelParam.Win_mean     = 162 *10^-12 ; % ~double Win_mean
        modelParam.del_G_sra    = 50* 3.0e-012;
        modelParam.inhSRA       = 0;
        modelParam.tau_sra = 100*10^(-3);	% SRA time constant (s)
        %}
        modelParam.vmInitMean  = -55.0e-3;
        disp('initialize_parameters(): PFtuned3 param set under revision')
    otherwise
        modelParam = override_params(modelParam, paramSet);
end


%% Set dependent parameters
modelParam = set_depedent_parameters(modelParam);


end