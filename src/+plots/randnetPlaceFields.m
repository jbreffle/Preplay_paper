function [] = randnetPlaceFields(modelParam, network, opS_PF, opV_PF, opGin_PF, linfields, varargin)
% plots.randnetPlaceFields(modelParam, network, opS_PF, opV_PF, opGin_PF, linfields, ithEnv=1, ithTrial=1)
%
% Plot results from place field simulations of randnet
%
% Inputs:
%   - modelParam: the model parameter structure
%   - network: the network structure
%   - opS_PF: logical spike matrix, opS(ithCell,ithTimeBin,ithTraj,ithTrial)
%   - opV_PF: membrane potential matrix, opV_PF(ithCell,ithTimeBin,ithTraj,ithTrial)
%   - opGin_PF: Input conductance matrix, opGin_PF(ithCell,ithTimeBin,ithTraj,ithTrial)
%   - linfields:
%
% Outputs:
%   - No outputs. Generates plots.
%
% Optional:
%   - ithEnv, ithTrial
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'network',	@isstruct)
addRequired(inputObj, 'opS_PF',     @islogical)
addRequired(inputObj, 'opV_PF',     @isnumeric)
addRequired(inputObj, 'opGin_PF', 	@isnumeric)
addRequired(inputObj, 'linfields', @iscell)
addParameter(inputObj, 'ithEnv',    1,  @isnumeric)
addParameter(inputObj, 'ithTrial',	1,	@isnumeric)
parse(inputObj, modelParam, network, opS_PF, opV_PF, opGin_PF, linfields, varargin{:});
p = inputObj.Results;

[~, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);

EtoPlot = network.E_indices(1:3);
ItoPlot = network.I_indices(1:3);


%% Histogram of mean firing rate over all trials of all environments
allMeanRates = sum(opS_PF, [2:4])./modelParam.t_max_PF./modelParam.nTrials_PF./modelParam.nEnvironments;
%allMeanRates_exampleEnv = sum(opS_PF(:,:,p.ithEnv,:), [2:4])./modelParam.t_max_PF./modelParam.nTrials_PF;
Erates = allMeanRates(network.E_indices);
Irates = allMeanRates(network.I_indices);
[~,BinEdges,~] = histcounts(allMeanRates, 50);

figure; hold on
xlabel('Mean rate (Hz, all PF trials)'); ylabel('Cell (count)')
histogram(Erates, BinEdges, FaceColor='r')
histogram(Irates, BinEdges, FaceColor='g')
legend({'E-cells', 'I-cells'})


%% Trial raster plots

myPlotSettings(width=8.5)
figure; tiledlayout(1,3, TileSpacing='compact')

%Example trial raster, unsorted
raster_ax(1) = nexttile;
title(['EnvID ', num2str(p.ithEnv), ', Trial ', num2str(p.ithTrial)], fontweight='normal')
plotSpikeRaster( logical( [ opS_PF(network.E_indices,:,p.ithEnv,p.ithTrial); opS_PF(network.I_indices,:,p.ithEnv,p.ithTrial) ]), 'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
xlabel('Time (s)'); ylabel('Cell');

% Same trial raster, E-cells sorted by PF peak location
raster_ax(2) = nexttile;
rpermIcells = randperm(numel(network.I_indices)); % Randomly permute cells, to prove they are not sorted
plotSpikeRaster( logical( [ opS_PF(network.E_indices(PFpeaksSequence{p.ithEnv}),:,p.ithEnv,p.ithTrial); opS_PF(network.I_indices(rpermIcells),:,p.ithEnv,p.ithTrial) ]), 'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
xlabel('Time (s)'); set(gca,'YTickLabel',[]); %ylabel('Cell');
title('Sorted by PF peak', fontweight='normal')

% Same trial raster, but sorted by just that trial
raster_ax(3) = nexttile;
rpermIcells = randperm(numel(network.I_indices)); % Randomly permute cells, to prove they are not sorted
trialESpikes = opS_PF(network.E_indices,:,p.ithEnv,p.ithTrial);
tmp = double(trialESpikes); tmp(tmp==0) = nan;
[~, trialESortInds] = sort(nanmean(tmp.*[1:size(tmp, 2)], 2), 'ascend');
plotSpikeRaster( logical( [trialESpikes(trialESortInds,:) ; opS_PF(network.I_indices(rpermIcells),:,p.ithEnv,p.ithTrial) ]), 'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
xlabel('Time (s)'); set(gca,'YTickLabel',[]); %ylabel('Cell');
title('Sorted by trial', fontweight='normal')

linkaxes(raster_ax, 'x')


%% Example V_m and G_in_PFs

myPlotSettings(lw=0.5)
figure; hold on

% E cell V_m
ex_ax(1) = nexttile;
plot(modelParam.t_PF, opV_PF(EtoPlot,:,p.ithEnv,p.ithTrial)*1e3);
box off
ylabel('V_m (mV)'); xlabel('Time (s)');
title('Example E-cells', fontweight='normal')

% I cell V_m
ex_ax(2) = nexttile;
plot(modelParam.t_PF, opV_PF(ItoPlot,:,p.ithEnv,p.ithTrial)*1e3);
box off
ylabel('V_m (mV)'); xlabel('Time (s)');
title('Example I-cells', fontweight='normal')

% E cell G_in
ex_ax(3) = nexttile;
plot(modelParam.t_PF, opGin_PF(EtoPlot,:,p.ithEnv,p.ithTrial)*1e9);
box off
ylabel(' G_{in} (pS)'); xlabel('Time (s)');
title('Example E-cell inputs', fontweight='normal')

% I cell G_in
ex_ax(4) = nexttile;
plot(modelParam.t_PF, opGin_PF(ItoPlot,:,p.ithEnv,p.ithTrial)*1e9);
box off
ylabel('G_{in} (pS)'); xlabel('Time (s)');
title('Example I-cell inputs', fontweight='normal')

linkaxes(ex_ax, 'x')
sgtitle(['EnvID ', num2str(p.ithEnv), ', Trial ', num2str(p.ithTrial)], fontSize=10)


%% Finish
drawnow
myPlotSettings


end