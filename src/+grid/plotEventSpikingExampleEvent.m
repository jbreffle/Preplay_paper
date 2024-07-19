function figHandle = plotEventSpikingExampleEvent(clusterMembershipMatrix, anyClusterIsActive, eventClustWiseRateSmoothed, eventRasterE, modelParam, analysisParam)
% figHandle = grid.plotEventSpikingExampleEvent(...);
%
% Plots example event with cluster-wise rates and raster plots.
% Called within
%
% TODO: Index the particular (paramSet, net, event) and directly plot that
% data (rather than calling within the grid.calculateEventSpiking()
% function)

%% Set up

if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=2, height=2.5)
end

MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
defaultLinesColormap = lines;


%% Plot event in different ways

% 1)
if analysisParam.plotExtraFigs
    [~, cellOrder] = sort(mean(clusterMembershipMatrix.*[1:modelParam.clusters]', 1));
    figure;
    tiledlayout(2,1, TileSpacing='tight', padding='tight');
    nexttile;
    hold on
    plot(eventClustWiseRateSmoothed);
    ylabel('Cluster rates (Hz)')
    h = gca;
    h.XAxis.Visible = 'off';
    box off
    nexttile;
    plotSpikeRaster( logical(eventRasterE(cellOrder,:)), ...
        'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
        'MarkerFormat', MarkerFormat);
    xlabel('Time (s)');
    ylabel('Cells (sorted)')
end

% 2) Only active cells, sorted only by active cluster order
if analysisParam.plotExtraFigs
    [~, mostActiveClust] = max(eventClustWiseRateSmoothed(anyClusterIsActive,:), [], 2);
    activeClusters = unique(mostActiveClust, 'stable');
    activeCells = any(eventRasterE, 2)';
    clusterOrder = nan(modelParam.clusters, 1); clusterOrder(activeClusters) = activeClusters;
    [~, cellOrder] = sort(nanmean(clusterMembershipMatrix.*clusterOrder, 1));
    activeCells = any(eventRasterE(cellOrder,:), 2)';
    cellOrder(~activeCells) = [];
    %colorVec = ones(numel(cellOrder), 3)';
    %MarkerFormat.Color = colorVec;
    figure;
    tiledlayout(2,1, TileSpacing='tight', padding='tight');
    nexttile;
    hold on
    plot(eventClustWiseRateSmoothed);
    ylabel('Cluster rates (Hz)')
    h = gca;
    h.XAxis.Visible = 'off';
    box off
    nexttile;
    plotSpikeRaster( logical(eventRasterE(cellOrder,:)), ...
        'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
        'MarkerFormat', MarkerFormat);
    xlabel('Time (s)');
    ylabel('Cells (sorted)')
end

% 3) Cluster rasters plotted seperately
if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=2.25, height=2.8)
end
event_tVec = 0:modelParam.dt:(modelParam.dt*(size(eventClustWiseRateSmoothed, 1)-1));
[~, mostActiveClust] = max(eventClustWiseRateSmoothed(anyClusterIsActive,:), [], 2);
activeClusters = unique(mostActiveClust, 'stable');
nTiles = 1+numel(activeClusters)+1;
figure;
tiledlayout(nTiles,1, TileSpacing='none', padding='tight');
nexttile;
hold on
plot(event_tVec, eventClustWiseRateSmoothed);
ylabel({'Cluster', 'rate (Hz)'})
h = gca;
h.XAxis.Visible = 'off';
box off
activeCells = any(eventRasterE, 2)';
xlim([0, max(event_tVec)])
%activeCells = ones(size(activeCells)); % Comment this line out to plot only active cells
for k=activeClusters'
    MarkerFormat.Color = defaultLinesColormap(k,:);
    clusterCellInds = clusterMembershipMatrix(k,:);
    nexttile;
    hold on
    plotSpikeRaster( logical(eventRasterE(clusterCellInds&activeCells,:)), ...
        'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
        'MarkerFormat', MarkerFormat);
    yline(0)
    ylabel({'Cluster', [num2str(find(k==activeClusters')), ' cells'] })
    h = gca;
    h.XAxis.Visible = 'off';
    box off
    set(gca,'YDir','normal')
    xlim([0, max(event_tVec)])
end
MarkerFormat.Color = 'k';
clusterCellInds = ~any(clusterMembershipMatrix(activeClusters,:), 1)&activeCells;
nexttile;
plotSpikeRaster( logical(eventRasterE(clusterCellInds,:)), ...
    'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
    'MarkerFormat', MarkerFormat);
h = gca;
xlabel('Time (s)');
ylabel({'Other', 'cells'})
set(gca,'YDir','normal')
xlim([0, max(event_tVec)])

% 4) Cluster rasters plotted seperately, with cells in multiple clusters plotted seperately
if analysisParam.plotExtraFigs
    if isequal(analysisParam.figSettings, 'manuscript')
        myPlotSettings(width=2.25, height=2.8)
    end
    event_tVec = 0:modelParam.dt:(modelParam.dt*(size(eventClustWiseRateSmoothed, 1)-1));
    [~, mostActiveClust] = max(eventClustWiseRateSmoothed(anyClusterIsActive,:), [], 2);
    activeClusters = unique(mostActiveClust, 'stable');
    nTiles = 1 + numel(activeClusters) + numel(activeClusters) + 1;
    figHandle = figure;
    tiledlayout(nTiles,1, TileSpacing='none', padding='tight');
    nexttile;
    hold on
    plot(event_tVec, eventClustWiseRateSmoothed);
    ylabel({'Cluster', 'rate (Hz)'})
    h = gca;
    h.XAxis.Visible = 'off';
    box off
    xlim([0, max(event_tVec)])
    includedCells = ones(1, size(eventRasterE, 1)); % Comment this line out to plot only active cells
    for ithActiveClust=1:numel(activeClusters)
        % Plot cells that are in both ithActiveClust and ithActiveClust-1
        if ithActiveClust>1
            currentClusterInds = activeClusters([ithActiveClust-1, ithActiveClust]);
            nonCurrentClusterInds = activeClusters; nonCurrentClusterInds([ithActiveClust-1, ithActiveClust])=[];
            MarkerFormat.Color = mean( defaultLinesColormap(currentClusterInds,:), 1);
            clusterCellInds = all(clusterMembershipMatrix(currentClusterInds,:), 1) & ~any(clusterMembershipMatrix(nonCurrentClusterInds,:), 1);
            nexttile;
            hold on
            plotSpikeRaster( logical(eventRasterE(clusterCellInds&includedCells,:)), ...
                'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
                'MarkerFormat', MarkerFormat);
            yline(0)
            %ylabel({'Cluster', [num2str(find(k==activeClusters')), ' cells'] })
            h = gca; h.XAxis.Visible = 'off'; box off
            set(gca,'YDir','normal')
            xlim([0, max(event_tVec)])
        end
        % Plot cells in ithActiveClust
        currentClusterInds = activeClusters(ithActiveClust);
        nonCurrentClusterInds = activeClusters(1:end ~= ithActiveClust);
        MarkerFormat.Color = mean( defaultLinesColormap(currentClusterInds,:), 1);
        clusterCellInds = clusterMembershipMatrix(currentClusterInds,:) & ~any(clusterMembershipMatrix(nonCurrentClusterInds,:), 1);
        nexttile;
        hold on
        plotSpikeRaster( logical(eventRasterE(clusterCellInds&includedCells,:)), ...
            'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
            'MarkerFormat', MarkerFormat);
        yline(0)
        %ylabel({'Cluster', [num2str(find(k==activeClusters')), ' cells'] })
        h = gca;
        h.XAxis.Visible = 'off';
        box off
        set(gca,'YDir','normal')
        xlim([0, max(event_tVec)])
    end
    ylabel({'Cells', '(by cluster membership)'})
    % Plot cells in multiple non-adjacent active clusters
    MarkerFormat.Color = mean( defaultLinesColormap(1:numel(activeClusters),:), 1);
    clusterCellInds = all(clusterMembershipMatrix(activeClusters,:), 1) | all(clusterMembershipMatrix([activeClusters(1),activeClusters(end)],:), 1);
    if numel(activeClusters)~=3
        warning('This plot only works on events with three active clusters, need to generalize')
    end
    nexttile;
    plotSpikeRaster( logical(eventRasterE(clusterCellInds,:)), ...
        'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
        'MarkerFormat', MarkerFormat);
    yline(0)
    %ylabel('Cells in additional multiple active clusters cells')
    h = gca;
    h.XAxis.Visible = 'off';
    box off
    h.YAxis.Visible = 'on';
    box off
    set(gca,'YDir','normal')
    xlim([0, max(event_tVec)])
    % Plot non-active cluster cells
    MarkerFormat.Color = 'k';
    clusterCellInds = ~any(clusterMembershipMatrix(activeClusters,:), 1)&includedCells;
    nexttile;
    plotSpikeRaster( logical(eventRasterE(clusterCellInds,:)), ...
        'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
        'MarkerFormat', MarkerFormat);
    h = gca;
    %ylabel('Other cells')
    xlabel('Time (s)');
    set(gca,'YDir','normal')
    xlim([0, max(event_tVec)])
end


%% Return plot settings to default

myPlotSettings

end