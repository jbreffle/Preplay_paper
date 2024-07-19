function [] = randnetPlaceFieldsManuscript(modelParam, network, opS_PF, opV_PF, opGin_PF, linfields, varargin)
% plots.randnetPlaceFields(modelParam, network, opS_PF, opV_PF, opGin_PF, linfields, ithEnv=1, ithTrial=1)
%
% Plot results from place field simulations of randnet for manuscript
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

% Add spikes
opV_PF(opV_PF==modelParam.V_reset) = 40e-3; % Add spikes up to 40 mV

%% Example trial raster, E-cells sorted by PF peak location

myPlotSettings(width=2, height=2)
figure;
rpermIcells = 1:numel(network.I_indices); %randperm(numel(network.I_indices)); % Randomly permute cells, to prove they are not sorted
plotSpikeRaster( logical( [ opS_PF(network.E_indices(PFpeaksSequence{p.ithEnv}),:,p.ithEnv,p.ithTrial); opS_PF(network.I_indices(rpermIcells),:,p.ithEnv,p.ithTrial) ]), 'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
xlabel('Time (s)'); ylabel('Cell (sorted)');


%% Example E cell V_m

opV_PF(opV_PF==modelParam.V_reset)=0.04;
yAxisRanges = [-70, -49];
%EtoPlot = network.E_indices(1:4);
%EtoPlot = network.E_indices(randi(numel(network.E_indices), 1, 4));
EtoPlot = [255, 156, 208, 465]; % 205

myPlotSettings(width=4.25, height=3)% (lw=0.5)
figure; hold on

ex_ax(1) = nexttile;
plot(modelParam.t_PF, opV_PF(EtoPlot(1),:,p.ithEnv,p.ithTrial)*1e3);
box off
ylabel('V_m (mV)');
xlabel('Time (s)');
ylim(yAxisRanges)
h = gca; h.XAxis.Visible = 'off';
%title('Example E-cells', fontweight='normal')

ex_ax(2) = nexttile;
plot(modelParam.t_PF, opV_PF(EtoPlot(2),:,p.ithEnv,p.ithTrial)*1e3);
box off
ylabel('V_m (mV)');
xlabel('Time (s)');
ylim(yAxisRanges)
h = gca; h.XAxis.Visible = 'off';
h = gca; h.YAxis.Visible = 'off';
%title('Example I-cells', fontweight='normal')

ex_ax(3) = nexttile;
plot(modelParam.t_PF, opV_PF(EtoPlot(3),:,p.ithEnv,p.ithTrial)*1e3);
box off
ylabel('V_m (mV)');
xlabel('Time (s)');
ylim(yAxisRanges)
%title('Example E-cell inputs', fontweight='normal')

ex_ax(4) = nexttile;
plot(modelParam.t_PF, opV_PF(EtoPlot(4),:,p.ithEnv,p.ithTrial)*1e3);
box off
ylabel('V_m (mV)');
xlabel('Time (s)');
ylim(yAxisRanges)
h = gca; h.YAxis.Visible = 'off';
%title('Example I-cell inputs', fontweight='normal')

linkaxes(ex_ax, 'x')
linkaxes(ex_ax, 'y')
ylim([-60, -40])
ylim([-55, -45])
% sgtitle(['EnvID ', num2str(p.ithEnv), ', Trial ', num2str(p.ithTrial)], fontSize=10)


%% Finish
drawnow
myPlotSettings


end