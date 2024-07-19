function spikeMat = randnet_calculator_memOpt_wGin(modelParam, network, ithEnv, isPFSim, varargin)
% randnet_calculator_memOpt_wGin.m
%
% Simulates network of recurrent LIF neurons. Should produce identical
% spikeMat as randnet_calculator.m but runs slightly faster and with much
% lower memory requirements due to not storing all conductance matrices.
%
% Input:
%   - parameters: A struct containing simulation parameters
%   - network: A struct containing network connectivity information
%
% Optional Input Parameters (parameter-value pairs):
%   - G_i: A matrix of size (parameters.n, parameters.t_steps) representing the external conductance input to each neuron over time. Default is zeros.
%   - V_m0: A column vector of size (parameters.n, 1) representing the initial membrane potential of each neuron. Default is zeros.
%   - G_sra0: A column vector of size (parameters.n, 1) representing the initial refractory conductance of each neuron. Default is zeros.
%
% Output:
%   - spikeMat: A logical matrix representing the spike activity of neurons over time
%
% Note: another option is to have a single simulation function with a
% flag inside the loop that can save the variable vectors at each time
% step. But this if statement slows the simulation by ~3%
%
% Adapated from randnet_calculator_memOpt.m,.
% This version does not take in a full G_in matrix. Instead, only the
% vector of G_in values at each time step is calculated in-place.

%% Parse varargin
inputObj = inputParser;
addRequired(inputObj, 'parameters', @isstruct)
addRequired(inputObj, 'network', @isstruct)
addRequired(inputObj, 'ithEnv', @isnumeric)
addRequired(inputObj, 'isPFSim', @islogical)
addParameter(inputObj, 'initialize', true,	@islogical);
addParameter(inputObj, 'V_m0', zeros(modelParam.n, 1), @isnumeric)
addParameter(inputObj, 'G_sra0', zeros(modelParam.n, 1), @isnumeric)
addParameter(inputObj, 'seed', -1,    @isnumeric);
parse(inputObj, modelParam, network, ithEnv, isPFSim, varargin{:})
p = inputObj.Results;

% Only set the rng if a seed was passed in
if p.seed~=-1
    rng(p.seed, 'twister');
end


%% Set up

if isPFSim
    simDuration = modelParam.t_max_PF;
else
    simDuration = modelParam.t_max_preplay;
end

nTimeSteps = (simDuration/modelParam.dt)+1;

% Pre-allocation, (single time-step, updated each increment)
G_syn_I_E = zeros(modelParam.n,1); % Conductance for pre-inhib to post-excit (S)
G_syn_E_E = zeros(modelParam.n,1); % Conductance for pre-excit to post-excit (S)
G_syn_I_I = zeros(modelParam.n,1); % Conductance for pre-inhib to post-inhib (S)
G_syn_E_I = zeros(modelParam.n,1); % Conductance for pre-excit to post-inhib (S)
G_sra = p.G_sra0; % Initial SRA vector
V_m = p.V_m0; % Initial Vm vector

% Initialize output spike matrix
spikeMat = false(modelParam.n, nTimeSteps);

% Binary indices of excitatory and inhibitory neurons
Einds = ismember(1:modelParam.n, network.E_indices)';
Iinds = ismember(1:modelParam.n, network.I_indices)';

% Pre-calculate constants:
expDecayE = exp(-modelParam.dt/modelParam.tau_syn_E);
expDecayI = exp(-modelParam.dt/modelParam.tau_syn_I);
expDecaySRA = exp(-modelParam.dt/modelParam.tau_sra);


%% G_in specific set-up

if isfield(modelParam, 'spatialInputType')
    inputType = modelParam.spatialInputType;
    if isequal(modelParam.spatialInputType, 'stepped')
        inputNSteps = modelParam.inputNSteps; % 5
    end
else % linear or stepped
    inputType = 'linear';
end

% Initialize G_in to approximate steady-state values (by default)
G_in = zeros(modelParam.n, 1);
if p.initialize
    G_in = 1/2 * modelParam.Win_mean * 2 * modelParam.rG * modelParam.tau_syn_E + ...
        sqrt(1/2*modelParam.tau_syn_E*modelParam.Win_mean.^2*2*modelParam.rG).*randn(modelParam.n, 1);
end

% Pre-compute:
if isPFSim % If place field simulation
    contextScalingE = Einds.*network.contextInput(:,ithEnv).* modelParam.PFcontextFrac;
    contextScalingI = Iinds.*network.contextInput(:,ithEnv).* modelParam.IcueScale_PF;
    spatialScaling1 = network.spatialInput{1}(:,ithEnv).*(1-modelParam.PFcontextFrac);
    spatialScaling2 = network.spatialInput{2}(:,ithEnv).*(1-modelParam.PFcontextFrac);
    if isequal(inputType, 'linearTernary')
        spatialScaling3 = network.spatialInput{3}(:,ithEnv).*(1-modelParam.PFcontextFrac);
    end
else % If preplay simulation
    if isfield(modelParam, 'useSleepContext') && modelParam.useSleepContext
        sleepContext = network.sleepContext;
    else
        sleepContext = network.contextInput(:,ithEnv);
    end
    contextScalingE = Einds.*sleepContext.* 1;
    contextScalingI = Iinds.*sleepContext.* modelParam.IcueScale ;
end

% Spatially modulated input functions:
dtBaseRate = modelParam.dt* modelParam.rG;
% Linearly scales across space, to peak at either end of track
linearInput = @(i) ...
    spatialScaling1.* (rand(modelParam.n, 1) < (dtBaseRate*(i/nTimeSteps))) + ...
    spatialScaling2.* (rand(modelParam.n, 1) < (dtBaseRate*((nTimeSteps-i)/nTimeSteps)));
% Three linearly varying inputs, third peaks at center
linearTernaryInput = @(i) ...
    spatialScaling1.* (rand(modelParam.n, 1) < (dtBaseRate* max(0, ((2*i)/(nTimeSteps)-1))              )) + ...
    spatialScaling2.* (rand(modelParam.n, 1) < (dtBaseRate* max(0, (2*(nTimeSteps-i)/nTimeSteps-1))     )) + ...
    spatialScaling3.* (rand(modelParam.n, 1) < (dtBaseRate* (1-(1/nTimeSteps*2)*abs(i-(nTimeSteps/2)))  ));
% Stepped input
steppedInput = @(i) spatialScaling1.* (rand(modelParam.n, 1)<(dtBaseRate*round(inputNSteps*i/nTimeSteps)/inputNSteps)) + ...
    spatialScaling2.* (rand(modelParam.n, 1)<(dtBaseRate*round(inputNSteps*(nTimeSteps-i)/nTimeSteps)/inputNSteps));

switch inputType
    case 'linear'
        spatialInputStep = linearInput;
    case 'linearTernary'
        spatialInputStep = linearTernaryInput;
    case 'stepped'
        spatialInputStep = steppedInput;
    otherwise
        error("create_G_in: unknown inputType option")
end

if ~isPFSim
    spatialInputStep = @(i) 0;
end

contextInputStep = @() contextScalingE .* (rand(modelParam.n, 1)<dtBaseRate) + ...
    contextScalingI .* (rand(modelParam.n, 1)<dtBaseRate);

%% Simulation
for tInd = 1:nTimeSteps % Need to go to nTimeSteps-1 if comparing to randnet_calculator.m

    % Check for spiking neurons
    spikes_all = [V_m>=modelParam.V_th];	% Boolean index of all cells spiking this timestemp
    spikers_I = [spikes_all & Iinds];       % Boolean vector of cells that spiked and are I cells
    spikers_E = [spikes_all & Einds];       % Boolean vector of cells that spiked and are E cells
    spikeMat(:,tInd) = spikes_all;          % Store spikes

    % Increment spiking SRA conductances
    if modelParam.inhSRA   % All cells have SRA
        G_sra = G_sra + spikes_all.*modelParam.del_G_sra;
    else        % Only E-cells have SRA
        G_sra = G_sra + spikers_E.*modelParam.del_G_sra;
    end

    % Increment synaptic conductances from spikes
    incoming_conn_E = sum(network.conns(spikers_E,:),1)'; % Post-synaptic neuron E input counts
    incoming_conn_I = sum(network.conns(spikers_I,:),1)'; % Post-synaptic neuron I input counts
    G_syn_I_E = G_syn_I_E + modelParam.del_G_syn_I_E*incoming_conn_I.*Einds;
    G_syn_E_E = G_syn_E_E + modelParam.del_G_syn_E_E*incoming_conn_E.*Einds;
    G_syn_I_I = G_syn_I_I + modelParam.del_G_syn_I_I*incoming_conn_I.*Iinds;
    G_syn_E_I = G_syn_E_I + modelParam.del_G_syn_E_I*incoming_conn_E.*Iinds;

    % Update membrane potential
    Gtot = modelParam.G_L + G_sra + G_syn_E_E + G_syn_E_I + G_syn_I_I + G_syn_I_E + G_in;
    V_ss = (G_in.*modelParam.V_syn_E + G_syn_E_E.*modelParam.V_syn_E + G_syn_E_I.*modelParam.V_syn_E + G_syn_I_I.*modelParam.V_syn_I + G_syn_I_E.*modelParam.V_syn_I + modelParam.G_L*modelParam.E_L + G_sra*modelParam.E_K)./Gtot;
    taueff = modelParam.C_m./Gtot;
    V_m = V_ss + (V_m - V_ss).*exp(-modelParam.dt ./taueff);

    % Update those that just spiked to reset
    V_m(spikes_all) = modelParam.V_reset;

    % Conductance exponential decays
    G_sra = G_sra.*expDecaySRA;      % SRA conductance update
    G_syn_E_E = G_syn_E_E.*expDecayE;	% E-E conductance update
    G_syn_I_E = G_syn_I_E.*expDecayI;	% I-E conductance update
    G_syn_I_I = G_syn_I_I.*expDecayI;	% I-I conductance update
    G_syn_E_I = G_syn_E_I.*expDecayE;	% E-I conductance update
    G_in = G_in*expDecayE; % Feedforward input conductance update

    % Add new feedforward input spikes for the next timestep
    G_in = G_in + spatialInputStep(tInd+1) + contextInputStep(); % Poisson spiking input functions

end
