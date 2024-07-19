%% Plot small example networks
% Clustered network and Watts-Strogatz networks with various beta values


%% Setup

addpath(fullfile("..", "..", "src"))
Config = utils.setConfig;

% Parameters and create network
analysisParam.figSettings = 'manuscript'; % standard, manuscript, SfNPoster
analysisParam.plotLegends = 1;
analysisParam.tSNESeed = 1234;

% Plot the network connectivity
switch analysisParam.figSettings
    case 'standard'
        myPlotSettings
    case 'manuscript'
        myPlotSettings(width=1.25, height=1)
    case 'SfNPoster'
        % myPlotSettings(width=3, height=1.5) % for poster
        myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3)
end


%% Small random clustered network

parameters.include_all = 2;
parameters.global_inhib = 1;

% Vary these parameters for figures
% Example manuscript figure is netSeed=1234, clusters=8, mnc=1.5
% Examples: [clusters, mnc] = [8, 1.5], [5, 1], [25, 2.0], [11, 3],
netSeed = 946; %randi(1000) % 1, 88, 946
analysisParam.tSNESeed = randi(100); % 1
parameters.clusters = 2;
parameters.mnc = 1.25;

parameters.n = 12;
parameters.p_I = 0.25;
parameters.n_I = parameters.n*parameters.p_I;

parameters.conn_prob = 0.2;

parameters.envIDs = 1;
parameters.Win_mu = -23.3568;
parameters.Win_sigma = 0.0694;
parameters.clusterCorrs = 1;
parameters.inputBiasSigma = 25;


parameters.p_E = 1-parameters.p_I;
parameters.npairs = parameters.n*(parameters.n-1);
parameters.cluster_n = round((parameters.mnc*parameters.n) / parameters.clusters) ;
parameters.nclusterpairs = parameters.cluster_n*(parameters.cluster_n - 1)*parameters.clusters;
parameters.cluster_prob = min(parameters.conn_prob*parameters.npairs/parameters.nclusterpairs,1);

% Create network
network = create_network(parameters, 'seed', netSeed);

rng(analysisParam.tSNESeed, 'twister')
% See plots.netDimRed.m for alternative plotting options
E_only = 1; % only plot E-E connections
dimRedInput = 'WW'; % 'W', 'WW', 'normClust', 'clust'
scatterSize = [1, 1];
scatterSize = 0.5*[1, 1];

network_indices = network.E_indices;
W = network.conns(network_indices, network_indices);
X = [W, W']; % to use both inputs and outputs
Y_tsne = tsne(X);

c1 = network.cluster_mat(2,network_indices)';
c2 = network.cluster_mat(1,network_indices)';
c = c1*[0, 1, 0] + c2*[1, 0, 0]; % color of scatter plot points
sz = scatterSize(1) + scatterSize(2)*sum(network.cluster_mat(:,network_indices))'; % size of scatter plot points
sz = sz*5;

%Wdigraph = digraph( W .* [any(c, 2)*any(c, 2)']);
Wdigraph = digraph( W); % .* [c1*c2']); % both pre and post synaptic neurons are in group 1
% Wdigraph = digraph( W .* [c1]); % presynaptic neuron is in group 1

% Plot t-SNE
figure
tmp = plot(Wdigraph, ... %, 'XData', Y_tsne(:,1),'YData', Y_tsne(:,2), ...
    'NodeLabel', [], 'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], ...
    'EdgeAlpha',0.5, 'LineWidth',1.0 ...
    );
% title(['tsne plot of [' dimRedInput, ']'])
axis off

randNetSWI = calcSWI(logical(W));


%% Plot empty figure with legend

if analysisParam.plotLegends
    figure; hold on; h = zeros(3, 1);
    h(1) = scatter(NaN,NaN,'r', 'filled');
    h(2) = scatter(NaN,NaN,'g', 'filled');
    h(3) = scatter(NaN,NaN,'y', 'filled');
    legend(h, 'Cluster 1', 'Cluster 2', 'Both clusters', 'Location', 'Best');
    axis off
end


%% Watts-Strogatz graph

rng(10)

% Plotting parameters
sz = 5;
c = 'k';

% Network parameters
N = 9;
p = 0.2;
beta = 0.00;

% No re-wiring
W = logical(full(adjacency(WattsStrogatz(N, round(N*p), beta))));
Wdigraph = digraph(W);
latticeSWI = calcSWI(logical(W));
figure
tmp = plot(Wdigraph, ...
    'NodeLabel', [], 'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], ...
    'EdgeAlpha',0.5, 'LineWidth',1.0 ...
    );
axis off

% To plot all on same x/y axes
xPos = tmp.XData;
yPos = tmp.YData;

% Partial re-wiring
beta = 0.2;
W = logical(full(adjacency(WattsStrogatz(N, round(N*p), beta))));
Wdigraph = digraph(W);
WattsStrogatzSWI = calcSWI(logical(W));
figure;
tmp = plot(Wdigraph, 'XData', xPos,'YData', yPos, ...
    'NodeLabel', [], 'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], ...
    'EdgeAlpha',0.5, 'LineWidth',1.0 ...
    );
axis off

% Complete re-wiring
beta = 1.0;
W = logical(full(adjacency(WattsStrogatz(N, round(N*p), beta))));
Wdigraph = digraph(W);
randomSWI = calcSWI(logical(W));
figure;
tmp = plot(Wdigraph, 'XData', xPos,'YData', yPos, ...
    'NodeLabel', [], 'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], ...
    'EdgeAlpha',0.5, 'LineWidth',1.0 ...
    );
axis off
