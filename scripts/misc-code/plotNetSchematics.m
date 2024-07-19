%% This code plots an example newtork figure
% plotNetSchematics.m
%

addpath(fullfile("..", "..", "src"))
Config = utils.setConfig;

% Parameters and create network
analysisParam.figSettings = 'manuscript'; % standard, manuscript, poster
analysisParam.plotLegends = 1;
analysisParam.tSNESeed = 1234;


%% Create network

% Vary these parameters for presentation figures:
% Example manuscript figure is netSeed=1234, clusters=8, mnc=1.5
% Examples for presentation: [clusters, mnc] = [8, 1.5], [5, 1], [25, 2.0], [11, 3],
netSeed = 1234;
parameters.clusters = 8;
parameters.mnc = 1.5;

parameters.n = 500;
parameters.n_I = 125;
parameters.p_I = 0.25;
parameters.conn_prob = 0.08;

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

parameters.include_all = 2;
parameters.global_inhib = 1;

% Create network
network = create_network(parameters, 'seed', netSeed);


%% Plot the network connectivity
switch analysisParam.figSettings
    case 'standard'
        myPlotSettings
    case 'manuscript'
        myPlotSettings(width=2, height=1.5)
    case 'poster'
        % myPlotSettings(width=3, height=1.5) % for poster
        myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3) % Poster format
end

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

%Wdigraph = digraph( W .* [any(c, 2)*any(c, 2)']);
Wdigraph = digraph( W .* [c1*c1']); % both pre and post synaptic neurons are in group 1
% Wdigraph = digraph( W .* [c1]); % presynaptic neuron is in group 1

% Plot t-SNE
figure
Xpos = Y_tsne(:,1); Ypos = Y_tsne(:,2);
plot(Wdigraph, 'XData', Xpos,'YData', Ypos, ...
    'MarkerSize',sz, 'NodeColor', c, 'EdgeColor', 0.3*[1 1 1], 'EdgeAlpha',0.5, 'LineWidth',1.0) % 'EdgeAlpha',0.5, 'LineWidth',1.0)
% title(['tsne plot of [' dimRedInput, ']'])
axis off

% Plot empty figure with legend
if analysisParam.plotLegends
    figure; hold on; h = zeros(4, 1);
    h(1) = scatter(NaN,NaN,'r', 'filled'); h(2) = scatter(NaN,NaN,'g', 'filled');
    h(3) = scatter(NaN,NaN,'y', 'filled'); h(4) = scatter(NaN,NaN,'k', 'filled');
    legend(h, ' 1, ~2', '~1,  2', ' 1,  2', '~1, ~2', 'Location', 'Best');
    axis off
end


%% Plot the network's histograms
switch analysisParam.figSettings
    case 'standard'
        myPlotSettings
    case 'manuscript'
        myPlotSettings(width=1.3, height=1.2)
    case 'poster'
        % myPlotSettings(width=3, height=1.5) % for poster
        myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3) % SfN-poster format
end

figure; histogram(sum(network.cluster_mat, 1)); xlabel('Clusters (count)'); ylabel('Cells (count)')

allInputs = vertcat(network.spatialInput{:});
figure; histogram(allInputs(allInputs~=0)*10^12); xlabel('Input strength (pS)'); ylabel('Inputs (count)')
% figure; histogram(sum(network.cluster_mat, 2)); xlabel('n Neurons'); ylabel('Clusters (count)')

EMat = network.conns(network.E_indices, network.E_indices);
figure; histogram(mean(EMat, 1), 30); xlabel('Frac. network'); ylabel('Inputs (count)')


%% Plot main input schematics
switch analysisParam.figSettings
    case 'standard'
        myPlotSettings
    case 'manuscript'
        myPlotSettings(width=2, height=1)
    case 'poster'
        % myPlotSettings(width=3, height=1.5) % for poster
        myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3) % SfN-poster format
end

% Plot Gaussians
xx = 0:0.1:10;
%yy = [gaussmf(xx,[1 2]); gaussmf(xx,[1 3]); gaussmf(xx,[1 4]); gaussmf(xx,[1 5]); gaussmf(xx,[1 6]); gaussmf(xx,[1 7]); gaussmf(xx,[1 8]);]
yy = [gaussmf(xx,[1 4]); gaussmf(xx,[1 5]); gaussmf(xx,[1 6]);];
figure; plot(xx,yy)
xlabel('Position'); ylabel('Input rate (Hz)')
box off; set(gca,'xtick',[]); set(gca,'ytick',[])
% legend({'Cell 3''s input', 'Cell 4''s input', 'Cell 5''s input'}); legend boxoff
% legend({'Input to cell 3', 'Input to cell 4', 'Input to cell 5'}); legend boxoff
ylim([min(yy(:)), 1.5*max(yy(:))])

% Plot linearly-varying cues
xx = 0:0.5:1;
yy = [ ([6, 6, 6].*xx(end)); (6.*xx); (6.*fliplr(xx))];
figure; plot(xx,yy')
xlabel('Position'); ylabel('Input rate (kHz)')
box off;
set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
ylim([min(yy(:)), 1.5*max(yy(:))])


%% Plot legends separately

if analysisParam.plotLegends
    figure; plot([0 0; 0 0; 0 0]'); legend({'Context cue', 'Location cue 1', 'Location cue 2'}, 'Box', 'off')
    axis off
    figure; plot([0 0; 0 0; 0 0]'); legend({'Input to cell i-1', 'Input to cell i', 'Input to cell i+1'}, 'Box', 'off')
    axis off
    figure; plot([0 0; 0 0; 0 0; 0 0]'); legend({'Context cue', 'Location cue 1', 'Location cue 2', 'Location cue 3'}, 'Box', 'off')
    axis off
end


%% Supplemental input schematics
switch analysisParam.figSettings
    case 'standard'
        myPlotSettings
    case 'manuscript'
        %myPlotSettings(width=3.0, height=2.0)
        myPlotSettings(width=1.0, height=0.8)
    case 'poster'
end

%figure

% Adding third cue that is center peaked
if 0
    xx = 0:0.5:1;
    yy = [ ([6, 6, 6].*xx(end)); (6.*xx); (6.*fliplr(xx)); ([0, 6, 0])];
    figure; plot(xx,yy')
    xlabel('Position'); ylabel('Input rate (kHz)')
    box off;
    set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
    ylim([min(yy(:)), 1.5*max(yy(:))])
end

% Adding third cue that is center peaked and equalizing total input
dt = 0.01;
tMax = 2;
t = 0:dt:tMax;
nTimeSteps = numel(t);
i = 1:nTimeSteps;
yy = [
    ones(size(i));
    max(0, ((2*i)/(nTimeSteps)-1));
    max(0, (2*(nTimeSteps-i)/nTimeSteps-1));
    1-(1/nTimeSteps*2)*abs(i-(nTimeSteps/2))
    ];
figure;
plot(t, yy)
xlabel('Position'); ylabel('Input rate (kHz)')
box off;
set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
ylim([min(yy(:)), 1.25*max(yy(:))])

% Plot stepped inputs
% 5 steps
xx=1:1:10000;
inputNSteps=5;
yy = [ones(size(xx)); round(inputNSteps*xx/numel(xx))/inputNSteps; round(inputNSteps*(numel(xx)-xx)/numel(xx))/inputNSteps];
figure
plot(xx, yy')
xlabel('Position'); ylabel('Input rate (kHz)')
box off;
set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
ylim([min(yy(:)), 1.25*max(yy(:))])

% 3 steps
xx=1:1:10000;
inputNSteps=3;
yy = [ones(size(xx)); round(inputNSteps*xx/numel(xx))/inputNSteps; round(inputNSteps*(numel(xx)-xx)/numel(xx))/inputNSteps];
figure
plot(xx, yy')
xlabel('Position'); ylabel('Input rate (kHz)')
box off;
set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
ylim([min(yy(:)), 1.25*max(yy(:))])

% 1 step
xx=1:1:10000;
inputNSteps=1;
yy = [ones(size(xx)); round(inputNSteps*xx/numel(xx))/inputNSteps; round(inputNSteps*(numel(xx)-xx)/numel(xx))/inputNSteps];
figure
plot(xx, yy')
xlabel('Position'); ylabel('Input rate (kHz)')
box off;
set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
ylim([min(yy(:)), 1.25*max(yy(:))])

% Plot linearly-varying cues
xx = 0:0.5:1;
yy = [ ([6, 6, 6].*xx(end)); (6.*xx); (6.*fliplr(xx))];
figure; plot(xx,yy')
xlabel('Position'); ylabel('Input rate (kHz)')
box off;
set(gca,'xtick',[]); set(gca,'ytick',[]); ylabel('Input rate (Hz)')
ylim([min(yy(:)), 1.25*max(yy(:))])
