function modelParam = set_depedent_parameters(modelParam, varargin)
% set_depedent_parameters.m
%
% Takes in "parameters" structure of initialized independent parameters and
% sets (or updates) all dependent parameters.
%
% Input:
%   A parameters structure
%
% Output:
%   The parameters structure with the dependent parameters created
%

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'modelParam',	@isstruct)
parse(inputObj, modelParam, varargin{:});
p = inputObj.Results;

%% Set the dependent parameters

modelParam.t_steps_PF = (modelParam.t_max_PF/modelParam.dt) + 1;	% Number of timesteps in simulation
modelParam.t_steps_preplay = (modelParam.t_max_preplay/modelParam.dt) + 1;	% Number of timesteps in simulation

% To make Wei and Wie equal, set one to nan before calling set_dependent_parameters
if isnan(modelParam.del_G_syn_I_E) && ~isnan(modelParam.del_G_syn_E_I)
    modelParam.del_G_syn_I_E = modelParam.del_G_syn_E_I;
elseif isnan(modelParam.del_G_syn_E_I) && ~isnan(modelParam.del_G_syn_I_E)
    modelParam.del_G_syn_E_I = modelParam.del_G_syn_I_E;
elseif isnan(modelParam.del_G_syn_E_I) && isnan(modelParam.del_G_syn_I_E)
    error('Wei and Wie are both nan')
end

% Calculate connection probabilites
modelParam.cluster_n = round((modelParam.mnc*modelParam.n)/modelParam.clusters);	% Average number of neurons in a cluster
modelParam.cluster_memb_prob = modelParam.mnc/modelParam.clusters;                  % Probability any given neuron is in any given cluster
modelParam.npairs = modelParam.n*(modelParam.n-1);                                  % Total number of possible neuron connections
modelParam.nclusterpairs = modelParam.cluster_n*(modelParam.cluster_n - 1)*modelParam.clusters;   % Total number of possible intra-cluster connections
modelParam.cluster_prob = min(modelParam.conn_prob*modelParam.npairs/modelParam.nclusterpairs,1); % Intra-cluster connection probability
modelParam.n_I = round((1-modelParam.p_E)*modelParam.n);	% Number of inhibitory neurons
modelParam.n_E = modelParam.n-modelParam.n_I;               % Number of excitatory neurons

% Set lognormal distribution mu and sigma, based on desired mean and variance
modelParam.Win_mu = log(modelParam.Win_mean^2 / sqrt(modelParam.Win_var+modelParam.Win_mean^2));
modelParam.Win_sigma = sqrt(log(modelParam.Win_var/modelParam.Win_mean^2 + 1));


% PF params
modelParam.nEnvironments = numel(modelParam.envIDs);
modelParam.gridxvals = modelParam.spatialBin:modelParam.spatialBin:modelParam.trackWidth; % Grid-points in the x-direction


end