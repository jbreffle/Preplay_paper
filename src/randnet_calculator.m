function [spikeMat, V_m, G_sra, G_syn_E_E, G_syn_I_E, G_syn_E_I, G_syn_I_I] = randnet_calculator(modelParam, network, simDuration, varargin)
% randnet_calculator.m
%
% Simulates network of recurrent LIF neurons
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
%   - V_m: A matrix representing the membrane potential of neurons over time
%   - G_sra: A matrix representing the refractory conductance of neurons over time
%   - G_syn_E_E: A matrix representing the conductance for pre-excitatory to post-excitatory synapses over time
%   - G_syn_I_E: A matrix representing the conductance for pre-inhibitory to post-excitatory synapses over time
%   - G_syn_E_I: A matrix representing the conductance for pre-excitatory to post-inhibitory synapses over time
%   - G_syn_I_I: A matrix representing the conductance for pre-inhibitory to post-inhibitory synapses over time

%% Parse varargin
inputObj = inputParser;
addRequired(inputObj, 'modelParam', @isstruct)
addRequired(inputObj, 'network', @isstruct)
addRequired(inputObj, 'simDuration', @isnumeric)
addParameter(inputObj, 'G_in', zeros(modelParam.n, (simDuration/modelParam.dt)+1), @isnumeric)
addParameter(inputObj, 'V_m0', zeros(modelParam.n, 1), @isnumeric)
addParameter(inputObj, 'G_sra0', zeros(modelParam.n, 1), @isnumeric)
parse(inputObj, modelParam, network, simDuration, varargin{:})
p = inputObj.Results;


%% Set up

simTSteps = (simDuration/modelParam.dt)+1;

% Pre-allocation
G_syn_I_E = zeros(modelParam.n,simTSteps);	% Conductance for pre-inhib to post-excit (S)
G_syn_E_E = zeros(modelParam.n,simTSteps);   % Conductance for pre-excit to post-excit (S)
G_syn_I_I = zeros(modelParam.n,simTSteps);   % Conductance for pre-inhib to post-inhib (S)
G_syn_E_I = zeros(modelParam.n,simTSteps);   % Conductance for pre-excit to post-inhib (S)
G_sra = zeros(modelParam.n,simTSteps);       % Refractory conductance for each neuron at each timestep (S)
G_sra(:,1) = p.G_sra0;  % Initial SRA vector
V_m = zeros(modelParam.n,simTSteps);
V_m(:,1) = p.V_m0;      % Initial V_m vector

spikeMat = false(modelParam.n, simTSteps);   % Initialize output spike matrix

% Binary indices of excitatory and inhibitory neurons
Einds = ismember(1:modelParam.n, network.E_indices)';
Iinds = ismember(1:modelParam.n, network.I_indices)';

% Pre-calculation:
expDecayE = exp(-modelParam.dt/modelParam.tau_syn_E);
expDecayI = exp(-modelParam.dt/modelParam.tau_syn_I);
expDecaySRA = exp(-modelParam.dt/modelParam.tau_sra);

%% Simulation
for tInd = 1:simTSteps-1

    % Check for spiking neurons
    spikes_all = [V_m(:,tInd) >= modelParam.V_th];	% Boolean index of all cells spiking this timestemp
    spikers_I = [spikes_all & Iinds];               % Boolean vector of cells that spiked and are I cells
    spikers_E = [spikes_all & Einds];               % Boolean vector of cells that spiked and are E cells
    spikeMat(:,tInd) = spikes_all;                  % Store spikes

    % Increment spiking SRA conductances
    if modelParam.inhSRA   % All cells have SRA
        G_sra(spikes_all,tInd) = G_sra(spikes_all,tInd) + modelParam.del_G_sra;
    else        % Only E-cells have SRA
        G_sra(spikers_E,tInd) = G_sra(spikers_E,tInd) + modelParam.del_G_sra;
    end

    % Increment synaptic conductances from spikes
    incoming_conn_E = sum(network.conns(spikers_E,:),1)'; % Post-synaptic neuron E input counts
    incoming_conn_I = sum(network.conns(spikers_I,:),1)'; % Post-synaptic neuron I input counts
    G_syn_I_E(:,tInd) = G_syn_I_E(:,tInd) + modelParam.del_G_syn_I_E*incoming_conn_I.*Einds;
    G_syn_E_E(:,tInd) = G_syn_E_E(:,tInd) + modelParam.del_G_syn_E_E*incoming_conn_E.*Einds;
    G_syn_I_I(:,tInd) = G_syn_I_I(:,tInd) + modelParam.del_G_syn_I_I*incoming_conn_I.*Iinds;
    G_syn_E_I(:,tInd) = G_syn_E_I(:,tInd) + modelParam.del_G_syn_E_I*incoming_conn_E.*Iinds;

    % Update membrane potential
    Gtot = (modelParam.G_L + G_sra(:,tInd) + G_syn_E_E(:,tInd) + G_syn_E_I(:,tInd) + G_syn_I_I(:,tInd) + G_syn_I_E(:,tInd) + p.G_in(:,tInd));
    V_ss = ( p.G_in(:,tInd).*modelParam.V_syn_E + G_syn_E_E(:,tInd).*modelParam.V_syn_E + G_syn_E_I(:,tInd).*modelParam.V_syn_E + G_syn_I_I(:,tInd).*modelParam.V_syn_I + G_syn_I_E(:,tInd).*modelParam.V_syn_I + modelParam.G_L*modelParam.E_L + G_sra(:,tInd)*modelParam.E_K )./Gtot;
    taueff = modelParam.C_m./Gtot;
    V_m(:,tInd+1) = V_ss + (V_m(:,tInd) - V_ss).*exp(-modelParam.dt ./taueff);

    % Update those that just spiked to reset
    V_m(spikes_all,tInd+1) = modelParam.V_reset;

    % Conductance exponential decays
    G_sra(:,tInd+1) = G_sra(:,tInd)*expDecaySRA;        % SRA conductance update
    G_syn_E_E(:,tInd+1) = G_syn_E_E(:,tInd).*expDecayE; % E-E conductance update
    G_syn_I_E(:,tInd+1) = G_syn_I_E(:,tInd).*expDecayI; % I-E conductance update
    G_syn_I_I(:,tInd+1) = G_syn_I_I(:,tInd).*expDecayI; % I-I conductance update
    G_syn_E_I(:,tInd+1) = G_syn_E_I(:,tInd).*expDecayE; % E-I conductance update
end