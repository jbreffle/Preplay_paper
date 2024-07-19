function [] = randnetPreplay(modelParam, network, V_m, G_in, network_spike_sequences, varargin)
% plots.randnetPreplay(modelParam, network, V_m, G_var(ithTest).G_in, network_spike_sequences)
%
% Plots basic simulation results from randnet.m
%
% Inputs:
%   - modelParam: the model parameter structure
%   - network: the network structure
%   - V_m: membrane potential matrix for a single trial
%   - G_in: Input conductance matrix for a single trial
%   - network_spike_sequence: results struct from detect_PBE
%
% Outputs:
%   - No output. Generates plots.
%
% Optional:
%   - ithTest: which preplay test to analyze (move to optional
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'network',	@isstruct)
addRequired(inputObj, 'V_m',	@ismatrix)
addRequired(inputObj, 'G_in',	@ismatrix)
addRequired(inputObj, 'network_spike_sequences',	@isstruct)
addParameter(inputObj, 'ithTest', 1,    @isnumeric);
addParameter(inputObj, 'sortMethod', 2,    @isnumeric); % 0=unsorted, 1=sort by largest event, 2=sort by mean over events:
parse(inputObj, modelParam, network, V_m, G_in, network_spike_sequences,  varargin{:});
p = inputObj.Results;


%% Set up

% Formatting options for plotSpikeRaster: uncomment one line below
MarkerFormat = struct; % use only this line to choose plotSpikeRaster defaults
%MarkerFormat.MarkerSize = 4; %MarkerFormat.Marker = '.';
%MarkerFormat.LineWidth = 1.5; MarkerFormat.Marker = '|';

% Find spikes
spikeMat = V_m>=modelParam.V_th;

% Time vector
t = [0:modelParam.dt:modelParam.t_max_preplay];


%% If there are any detected events, plot sorted rasters

if isfield(network_spike_sequences(p.ithTest), 'events') && ...
        ~isempty(network_spike_sequences(p.ithTest).events)

    % Plot all sorted events
    myPlotSettings(width=8.5)
    figure
    events = network_spike_sequences(p.ithTest).events;
    num_events = size(events, 1);
    for ithEvent = 1:num_events

        spike_ranks = network_spike_sequences(p.ithTest).ranks_vec(:,ithEvent);
        if numel(spike_ranks)== modelParam.n_E % ranks for only E cells
            [~, Ie] = sort(spike_ranks);
            eventSpikes = [spikeMat(network.E_indices(Ie),events(ithEvent,1):events(ithEvent,2)); ...
                spikeMat(network.I_indices,events(ithEvent,1):events(ithEvent,2))];
        elseif numel(spike_ranks)==modelParam.n              % ranks for all cells
            [~, Ie] = sort(spike_ranks(network.E_indices));
            [~, Ii] = sort(spike_ranks(network.I_indices));
            eventSpikes = [spikeMat(network.E_indices(Ie),events(ithEvent,1):events(ithEvent,2)); ...
                spikeMat(network.I_indices(Ii),events(ithEvent,1):events(ithEvent,2))];
        else
            error('Unkonwn number of cells in spike_ranks')
        end
        subplot(1,num_events,ithEvent)
        plotSpikeRaster( eventSpikes, 'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat);
        title(num2str(ithEvent),'FontWeight','Normal')
        xlabel('Time (s)')
        if ithEvent == 1
            ylabel('Neuron (sorted)')
        else
            set(gca,'ytick',[])
        end

    end
    sgtitle('All sorted events')
    myPlotSettings
end


%% Plot example Gin and Vm

EtoPlot = network.E_indices(1:3);
ItoPlot = network.I_indices(1:3);

myPlotSettings(lw=0.5, width=8.5)
figure; hold on

% E cell V_m
ex_ax(1) = nexttile;
plot(t, V_m(EtoPlot,:)*1e3);
box off
ylabel('V_m (mV)'); xlabel('Time (s)');
title('Example E-cells', fontweight='normal')

% I cell V_m
ex_ax(2) = nexttile;
plot(t, V_m(ItoPlot,:)*1e3);
box off
ylabel('V_m (mV)'); xlabel('Time (s)');
title('Example I-cells', fontweight='normal')

% E cell G_in
ex_ax(3) = nexttile;
plot(t, G_in(EtoPlot,:)*1e9);
box off
ylabel(' G_{in} (pS)'); xlabel('Time (s)');
title('Example E-cell inputs', fontweight='normal')

% I cell G_in
ex_ax(4) = nexttile;
plot(t, G_in(ItoPlot,:)*1e9);
box off
ylabel('G_{in} (pS)'); xlabel('Time (s)');
title('Example I-cell inputs', fontweight='normal')

linkaxes(ex_ax, 'x')


%% All-cell raster with E and I pop. rates

myPlotSettings(width=8.5, height=5.5)

figure;
ax1 = subplot(2,1,1); hold on
if exist('events', 'var')
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, size(spikeMat, 1), size(spikeMat, 1)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
plotSpikeRaster( [spikeMat(network.E_indices,:);  spikeMat(network.I_indices,:)], ...
    'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat); %
ylabel('Cell'); xlabel('Time (s)');

% Plot E and I cell population firing rate
meanEPopRate = smoothdata(mean(spikeMat(network.E_indices,:), 1)/modelParam.dt, 'gaussian', modelParam.PBE_window);
meanIPopRate = smoothdata(mean(spikeMat(network.I_indices,:), 1)/modelParam.dt, 'gaussian', modelParam.PBE_window);
PBEthresh = max(mean(meanEPopRate)+(modelParam.PBE_zscore*std(meanEPopRate)), modelParam.PBE_min_Hz); % threshold rate for PBE detection

ax2 = subplot(2,1,2); hold on
if exist('events', 'var')
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, ceil(max(meanEPopRate)), ceil(max(meanEPopRate))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end

yyaxis right;
plot(t, meanIPopRate, ':')
ylabel('I cell mean rate (Hz)');
yyaxis left;
plot(t, meanEPopRate)
ylabel('E cell mean rate (Hz)'); xlabel('Time (s)');

yline(mean(meanEPopRate), 'g')
yline(mean(meanEPopRate)+ std(meanEPopRate))
yline(mean(meanEPopRate)+ 2*std(meanEPopRate))
yline(PBEthresh, 'r');

linkaxes([ax1, ax2], 'x')


%% E-cell raster and E-cell pop. rate
myPlotSettings(width=8.5, height=5.5)

figure;
if exist('events', 'var')
    sgtitle({['Frac. particip.: ', ...
        regexprep(num2str( round( mean(~isnan(network_spike_sequences(p.ithTest).ranks_vec), 1),  2)),'\s+',',') ], ...
        ['Event dur. (ms): ', ...
        regexprep(num2str( round( network_spike_sequences(p.ithTest).event_lengths'*1000,  0)),'\s+',', ') ] }, ...
        FontWeight='Normal', fontSize=12, HorizontalAlignment='center')
    [~, largestEventInd] = max(sum(~isnan(network_spike_sequences(p.ithTest).ranks_vec), 1));
    spike_ranks1 = network_spike_sequences(p.ithTest).ranks_vec(:,largestEventInd);
    [~, Ie1] = sort(spike_ranks1);
    meanRelRank = nanmean(network_spike_sequences(p.ithTest).ranks_vec./modelParam.n_E, 2);
    [~, IeRel] = sort(meanRelRank);
else
    sgtitle('No events detected')
    Ie1 = ones(size(network.E_indices));
end

if p.sortMethod==1 && exist('events', 'var')
    reordered_spikes = [spikeMat(network.E_indices(Ie1),:)];
elseif p.sortMethod==2 && exist('events', 'var')
    reordered_spikes = [spikeMat(network.E_indices(IeRel),:)];
else
    reordered_spikes = [spikeMat(network.E_indices,:)];
end

% Plot Raster
ax1 = subplot(2,1,1); hold on
if exist('events', 'var')
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, size(spikeMat, 1), size(spikeMat, 1)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
plotSpikeRaster( reordered_spikes, ...
    'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat); %
ylabel('E Cell (sorted)');

% Plot population firing rate
ax2 = subplot(2,1,2); hold on
meanEPopRate = smoothdata(mean(spikeMat(network.E_indices,:), 1)/modelParam.dt, 'gaussian', modelParam.PBE_window);
PBEthresh = max(mean(meanEPopRate)+(modelParam.PBE_zscore*std(meanEPopRate)), modelParam.PBE_min_Hz); % threshold rate for PBE detection
if exist('events', 'var')
    for i = 1:size(events, 1)
        fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, ceil(max(meanEPopRate)), ceil(max(meanEPopRate))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    end
end
ylim([0,  ceil(max(meanEPopRate))])
plot(t, meanEPopRate)
ylabel('E cell mean rate (Hz)'); xlabel('Time (s)');

yline(mean(meanEPopRate), '-g', {'+0z'})
yline(mean(meanEPopRate)+ std(meanEPopRate), '-k', {'+1z'})
yline(mean(meanEPopRate)+ 2*std(meanEPopRate), '-k', {'+2z'})
yline(PBEthresh, '--r', {'Threshold'}, LabelHorizontalAlignment='left');
linkaxes([ax1, ax2], 'x')
xlim([0,  t(end)*1.025])


%% Finish
drawnow
myPlotSettings


end