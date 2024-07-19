function clusterRates(spikeMat, modelParam, network, varargin)
% plots.clusterRates(spikeMat, modelParam, network);
%
% Generates a plot of the mean spike rate, mean cluster-wise spike rate,
% and z-scored cluster-wise spike rate from the data in spikes_V_m
%
% Inputs:
%	- spikeMat: E-cell binary spike matrix from one trial
%	- modelParam: the parameter structure for the simulation
%	- network: network structure
%
% % Outputs:
%
%
% Example usage, after running a simulation from randnet.m:
% smoothWindow = 10 * (1/parameters.dt * 1/1000); %gaussian kernel width for smoothing firing rate curves
% plots.clusterRates(spikes_V_m, parameters, network, 'smoothWindow', smoothWindow);

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'spikeMat',	@ismatrix)
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'network',	@isstruct)
addParameter(inputObj, 'smoothWindow', 1000,    @isnumeric);
parse(inputObj, spikeMat, modelParam, network, varargin{:});
p = inputObj.Results;

%% Default parameters:
legend_flag = 0; % if 1, include cluster legend
smoothWindow = 50 * (1/modelParam.dt * 1/1000); %gaussian kernel width for smoothing firing rate curves


%% Read in optional parameters, to overwrite above defaults
for i=1:2:length(varargin)
    switch varargin{i}
        case 'legend_flag'
            legend_flag = varargin{i+1};
        case 'smoothWindow'
            smoothWindow = varargin{i+1};
        otherwise
            error('plots.clusterRates: Unknown input')
    end
end


%% Main:
assert([size(spikeMat, 1)==modelParam.n_E], 'spikes_V_m must contain all E-cells')
t = [0:modelParam.dt:modelParam.t_max_preplay];

y = zeros(modelParam.clusters, size(spikeMat, 2)); % num spikes each cluster fired each time step
for ithCluster = 1:modelParam.clusters
    clusterMember = network.cluster_mat(ithCluster,network.E_indices);
    y(ithCluster,:) = clusterMember*spikeMat/sum(clusterMember)/modelParam.dt;
end

figure;

yRate = smoothdata(mean(spikeMat, 1)/modelParam.dt, 'gaussian', smoothWindow);
ax1 = subplot(3, 1, 1); plot(t, yRate);
ylabel('Pop. mean rate (Hz)')

ySmoothed = smoothdata(y, 2, 'gaussian', smoothWindow);
ax2 = subplot(3, 1, 2); plot(t, ySmoothed)
ylabel('Cluster mean rate (Hz)')
if legend_flag
    legend( sprintfc('Cluster %g', 1:modelParam.clusters), 'Location', 'Best' )
end

yZScore = (ySmoothed-mean(ySmoothed, 2))./std(ySmoothed, [], 2) ; % z-score across time, for each cluster
% yZScore = (ySmoothed-mean(ySmoothed, 1))./std(ySmoothed, [], 1) ; % z-score across clusters, for each time point
ax3 = subplot(3, 1, 3); plot(t, yZScore);;
ylabel('Cluster z-score')
xlabel('Time (s)');
if legend_flag
    legend( sprintfc('Cluster %g', 1:modelParam.clusters), 'Location', 'Best' )
end

linkaxes([ax1 ax2 ax3],'x')

end