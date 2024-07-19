function [] = updateEventSpikingPlots(paramGridFigHandles, analysisParam, gridOP, paramOP, simParam)
%UPDATEEVENTSPIKINGPLOTS Summary of this function goes here
%
% Updates the plots in the parameter grid loop for the event Spiking analysis script


ithGridOPPlot = 1;

xVals = simParam.variedParam(1).range;
xLabel = simParam.variedParam(1).name;
yVals = simParam.variedParam(2).range;
yLabel = simParam.variedParam(2).name;


% Spikes/cell/event analyses
if analysisParam.runSpikeDistAnalyses
    set(0,'CurrentFigure', paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.spikeDistMean', 'alphadata', ~isnan(gridOP.spikeDistMean')); set(gca,'YDir','normal'); colorbar; title('Mean spikes/cell/PBE')
    xlabel(xLabel); ylabel(yLabel)

    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.spikeDistCV', 'alphadata', ~isnan(gridOP.spikeDistMean'));
    set(gca,'YDir','normal');
    colorbar;
    title('C.V. of spikes/cell/PBE')
    xlabel(xLabel); ylabel(yLabel)
    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.spikeDistFitR2', 'alphadata', ~isnan(gridOP.spikeDistFitR2'));
    set(gca,'YDir','normal');
    colorbar;
    title('Exp fit R^2')
    xlabel(xLabel); ylabel(yLabel)
end

% Cluster-wise analyses
if analysisParam.runAllClustAnalyses
    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.clustSpikeRMSD', 'alphadata', ~isnan(gridOP.clustSpikeRMSD'));
    set(gca,'YDir','normal');
    colorbar;
    title('Mean SpikeClustRMSD')
    xlabel(xLabel);
    ylabel(yLabel)

    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.clustActivRMSD', 'alphadata', ~isnan(gridOP.clustActivRMSD'));
    set(gca,'YDir','normal');
    colorbar;
    title('Mean ClustActivRMSD')
    xlabel(xLabel);
    ylabel(yLabel)

    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.nClustMostActive', 'alphadata', ~isnan(gridOP.nClustMostActive')); set(gca,'YDir','normal'); colorbar; title('nClustActive')
    xlabel(xLabel);
    ylabel(yLabel)

    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.propMostActive', 'alphadata', ~isnan(gridOP.propMostActive'));
    set(gca,'YDir','normal');
    colorbar;
    title('Prop MostActive')
    xlabel(xLabel);
    ylabel(yLabel)

    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.clustMostActiveDur', 'alphadata', ~isnan(gridOP.clustMostActiveDur'));
    set(gca,'YDir','normal');
    colorbar;
    title('clustMostActiveDur')
    xlabel(xLabel);
    ylabel(yLabel)

    % Duration of active clusters
    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.clustActiveDur', 'alphadata', ~isnan(gridOP.clustActiveDur'));
    set(gca,'YDir','normal');
    colorbar;
    title('clustActiveLengths')
    xlabel(xLabel);
    ylabel(yLabel)

    % Most-active analyses
    set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
    ithGridOPPlot = ithGridOPPlot+1;
    imagesc(xVals, yVals, gridOP.fracMostActive', 'alphadata', ~isnan(gridOP.fracMostActive'));
    set(gca,'YDir','normal');
    colorbar; title('Frac. activity from most active')
    xlabel(xLabel);
    ylabel(yLabel)

    if size(gridOP.fracMostActive, 2) == numel(yVals)
        foldIncreaseFracMostActive = gridOP.fracMostActive ./ (1./yVals);
        set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
        ithGridOPPlot = ithGridOPPlot+1;
        imagesc(xVals, yVals, foldIncreaseFracMostActive', 'alphadata', ~isnan(foldIncreaseFracMostActive'));
        set(gca,'YDir','normal');
        colorbar;
        title('Fold. Frac. activity from most active')
        xlabel(xLabel);
        ylabel(yLabel)
    end
end

% Active cluster analyses:
% Expanded cluster activation window
set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
ithGridOPPlot = ithGridOPPlot+1;
imagesc(xVals, yVals, gridOP.clustActiveDur_expanded', 'alphadata', ~isnan(gridOP.clustActiveDur_expanded'));
set(gca,'YDir','normal');
colorbar;
title('Active cluster duration (s)','FontWeight','Normal')
xlabel(xLabel);
ylabel(yLabel)

% Number of clusters active
set(0,'CurrentFigure',paramGridFigHandles(ithGridOPPlot));
ithGridOPPlot = ithGridOPPlot+1;
imagesc(xVals, yVals, gridOP.nClustActive', 'alphadata', ~isnan(gridOP.nClustActive'));
set(gca,'YDir','normal');
colorbar;
title('Active clusters','FontWeight','Normal')
xlabel(xLabel);
ylabel(yLabel)


drawnow

end

