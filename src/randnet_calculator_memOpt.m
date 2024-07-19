function spikeMat = randnet_calculator_memOpt(modelParam, network, simDuration, varargin)
% randnet_calculator_memOpt.m
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


%% Parse varargin
inputObj = inputParser;
addRequired(inputObj, 'parameters', @isstruct)
addRequired(inputObj, 'network', @isstruct)
addRequired(inputObj, 'simDuration', @isnumeric)
addParameter(inputObj, 'G_in', zeros(modelParam.n, (simDuration/modelParam.dt)+1), @isnumeric)
addParameter(inputObj, 'V_m0', zeros(modelParam.n, 1), @isnumeric)
addParameter(inputObj, 'G_sra0', zeros(modelParam.n, 1), @isnumeric)
parse(inputObj, modelParam, network, simDuration, varargin{:})
p = inputObj.Results;


%% Set up

nTimeSteps = (simDuration/modelParam.dt)+1;

% Pre-allocation, (single time-step, updated each increment)
G_syn_I_E = zeros(modelParam.n,1);  % Conductance for pre-inhib to post-excit (S)
G_syn_E_E = zeros(modelParam.n,1);  % Conductance for pre-excit to post-excit (S)
G_syn_I_I = zeros(modelParam.n,1);  % Conductance for pre-inhib to post-inhib (S)
G_syn_E_I = zeros(modelParam.n,1);  % Conductance for pre-excit to post-inhib (S)
G_sra = p.G_sra0;   % Initial SRA vector
V_m = p.V_m0;       % Initial Vm vector

spikeMat = false(modelParam.n, nTimeSteps);   % Initialize output spike matrix

% Binary indices of excitatory and inhibitory neurons
Einds = ismember(1:modelParam.n, network.E_indices)';
Iinds = ismember(1:modelParam.n, network.I_indices)';

% Pre-calculation:
expDecayE = exp(-modelParam.dt/modelParam.tau_syn_E);
expDecayI = exp(-modelParam.dt/modelParam.tau_syn_I);
expDecaySRA = exp(-modelParam.dt/modelParam.tau_sra);

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
    Gtot = modelParam.G_L + G_sra + G_syn_E_E + G_syn_E_I + G_syn_I_I + G_syn_I_E + p.G_in(:,tInd);
    V_ss = (p.G_in(:,tInd).*modelParam.V_syn_E + G_syn_E_E.*modelParam.V_syn_E + G_syn_E_I.*modelParam.V_syn_E + G_syn_I_I.*modelParam.V_syn_I + G_syn_I_E.*modelParam.V_syn_I + modelParam.G_L*modelParam.E_L + G_sra*modelParam.E_K)./Gtot;
    taueff = modelParam.C_m./Gtot;
    V_m = V_ss + (V_m - V_ss).*exp(-modelParam.dt ./taueff);

    % Update those that just spiked to reset
    V_m(spikes_all) = modelParam.V_reset;

    % Conductance exponential decays
    G_sra =     G_sra.*expDecaySRA;      % SRA conductance update
    G_syn_E_E = G_syn_E_E.*expDecayE;	% E-E conductance update
    G_syn_I_E = G_syn_I_E.*expDecayI;	% I-E conductance update
    G_syn_I_I = G_syn_I_I.*expDecayI;	% I-I conductance update
    G_syn_E_I = G_syn_E_I.*expDecayE;	% E-I conductance update
end