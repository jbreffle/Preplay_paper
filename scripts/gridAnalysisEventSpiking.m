%% Analyze and plot within-event spiking and decode statisics across parameter grids
% gridAnalysisEventSpiking.m
%

%% Choose which simulation results to analyze

decodeFileID = "mainGrid300s"; % Select which decoding file result ID to analyze
%decodeFileID = "mainPaperGrid"; % Select which decoding file result ID to analyze


%% Set up and load data

myPlotSettings

addpath(fullfile("..", "src"))
Config = utils.setConfig;
dataPath = Config.simDataPath;

% Load parameters and results
decodeName = grid.getDecodeFilename(decodeFileID);
[modelParam, ...
    simParam, ...
    resultsStruct, ...
    PFresultsStruct, ...
    preplayGridSpikes, ...
    pfGridspikes ...
    ] = grid.loadResults(decodeName);


%% Analysis parameters
analysisParam.paramPoint = [2, 3]; % Select which example param point to plot, can be "all"
analysisParam.analyzeEntireGrid = false; % If false, only analyze analysisParam.paramPoint
analysisParam.ithEnv = 1; % which environment's decoding and place fields to plot/analyze
analysisParam.exampleEvent = [2, 3, 1, 1]; % Indices of example event to plot (param1, param2, net, event)

if analysisParam.analyzeEntireGrid
    analysisParam.paramSetInds = combvec([1:size(resultsStruct, 1)], [1:size(resultsStruct, 2)])';
else
    analysisParam.paramSetInds = analysisParam.paramPoint;
end

% Analysis parameters
analysisParam.minDurClustOn = 5e-3;
analysisParam.nEventRateBins = 100;
analysisParam.normEventLength = true;
analysisParam.smoothWindow = round(10 * 1e-3 / modelParam.dt); %round(1/10 * numel(event_timeVec));
analysisParam.activeClustMuliplyer = 2; % Clusters must be this times more active than any other cluster to be called 'active'
analysisParam.nShuffles = 100;

% Options to run/plot
analysisParam.plotExtraFigs = false;
analysisParam.runRasterChecks = false; % For each event, verify that the participating cells in the raster match the decode cells
analysisParam.runSpikeDistAnalyses = false; % Run analyses for the spikes/cell/PBE histograms
analysisParam.runAllClustAnalyses = false; % If false, only runs the 'active' cluster analyses
analysisParam.runShuffleFitting = false; % Analyses fitting linear regression to the shuffled events

analysisParam.figSettings = 'manuscript'; % 'manuscript' or ''

disp(analysisParam)
%assert(any(strcmp(analysisParam.figSettings, {'standard', 'manuscript', 'SfNPoster', 'manuscriptAlt'})))


%% Set up for analysis and plotting

gridOP = struct;
nGridFigs = 2; % Base number for manuscript only figures

% Spikes/cell/event analyses
if analysisParam.runSpikeDistAnalyses
    gridOP.spikeDistMean = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    gridOP.spikeDistCV = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    gridOP.spikeDistFitR2 = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    nGridFigs = nGridFigs + 3;
end

% Cluster-wise analyses
if analysisParam.runAllClustAnalyses
    gridOP.clustSpikeRMSD = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    gridOP.clustActivRMSD = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    gridOP.nClustMostActive = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    gridOP.propMostActive = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    gridOP.clustMostActiveDur = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2)));
    gridOP.clustActiveDur = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2))); % Mean duration that a cluster is active
    gridOP.fracMostActive = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2))); % Fraction of total activity due to the most active cluster
    nGridFigs = nGridFigs + 8;
end

% Main manuscript analyses
gridOP.clustActiveDur_expanded = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2))); % Mean duration that a cluster is active
gridOP.nClustActive = nan(max(analysisParam.paramSetInds(:,1)), max(analysisParam.paramSetInds(:,2))); % Number of clusters active per event

% Generate all figures now, so they can be updated in the background
if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=2.25, height=1.25)
end
paramGridFigHandles = nan(1, numel(fieldnames(gridOP)));
for i = 1:nGridFigs
    paramGridFigHandles(i) = figure;
end


if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=2.5, height=2.0)
else
    myPlotSettings
end


%% Parameter grid loop: analysis and plotting

tic
rng('default');
for ithParamSetToPlot = 1:size(analysisParam.paramSetInds, 1)

    % Set up varied-parameters for this loop
    param1Ind = analysisParam.paramSetInds(ithParamSetToPlot,1);
    param1Val = simParam.variedParam(1).range(param1Ind);
    disp(['Param1=', num2str(param1Ind), ': ', simParam.variedParam(1).name, '=', num2str(param1Val)])
    param2Ind = analysisParam.paramSetInds(ithParamSetToPlot,2);
    param2Val = simParam.variedParam(2).range(param2Ind);
    disp(['Param2=', num2str(param2Ind), ': ', simParam.variedParam(2).name, '=', num2str(param2Val)])
    linearParamInd = find(all([param1Val; param2Val] == simParam.parameterSets_vec, 1));
    if isempty(linearParamInd)
        continue
    end

    % Update parameters for the ithParamSetToPlot:
    loopModelParam = modelParam;
    loopModelParam.(simParam.variedParam(1).name) = simParam.variedParam(1).range(param1Ind);
    loopModelParam.(simParam.variedParam(2).name) = simParam.variedParam(2).range(param2Ind);
    loopModelParam = set_depedent_parameters(loopModelParam);

    % Run analysis calculation
    paramOP = grid.calculateEventSpiking(...
        analysisParam, simParam, loopModelParam, resultsStruct, PFresultsStruct, preplayGridSpikes{linearParamInd}, param1Ind, param2Ind...
        );

    % Store results
    if analysisParam.runSpikeDistAnalyses
        gridOP.spikeDistMean(param1Ind, param2Ind) = mean([paramOP.spikesPerPBE{:}], 'all');
        gridOP.spikeDistCV(param1Ind, param2Ind) =  std([paramOP.spikesPerPBE{:}], [], 'all') / mean([paramOP.spikesPerPBE{:}], 'all');
        gridOP.spikeDistFitR2(param1Ind, param2Ind) = paramOP.histFitR2;
    end
    if analysisParam.runAllClustAnalyses
        gridOP.clustSpikeRMSD(param1Ind, param2Ind) = mean(paramOP.clustSpikeRMSD);
        gridOP.clustActivRMSD(param1Ind, param2Ind) = mean([paramOP.eventClustRMSD{:}], 'all');
        gridOP.nClustMostActive(param1Ind, param2Ind) = mean(paramOP.nClustMostActive);
        gridOP.propMostActive(param1Ind, param2Ind) = nanmean([paramOP.event_prop_MostActClust{:}], 'all');
        gridOP.clustMostActiveDur(param1Ind, param2Ind) = mean(paramOP.clustMostActiveDur);
        gridOP.clustActiveDur(param1Ind, param2Ind) = mean(paramOP.clustActiveLengths);
        gridOP.fracMostActive(param1Ind, param2Ind) = nanmean(paramOP.fracMostActive);
    end
    gridOP.clustActiveDur_expanded(param1Ind, param2Ind) = mean(paramOP.clustActiveLengths_expanded);
    gridOP.nClustActive(param1Ind, param2Ind) = mean(paramOP.nClustActive);

    % Cluster co-activity anaylsis
    gridOP.coactiveSymmetryIndex(param1Ind, param2Ind) = mean(paramOP.coactiveSymmetryIndex);
    gridOP.coactiveSymmetryPVal(param1Ind, param2Ind) = median(paramOP.coactiveSymmetryPVal);

    % Correlation between cluster activation and decode-bin quality
    x = paramOP.binRMSE(:,4); y = paramOP.binRMSE(:,2); ind = ~isnan(x);
    mdl1 = fitlm(x(ind), y(ind));
    gridOP.mdl1PVal(param1Ind, param2Ind) = mdl1.Coefficients.pValue("x1");
    gridOP.mdl1Slope(param1Ind, param2Ind) = mdl1.Coefficients.Estimate("x1");
    x = paramOP.binRMSE(:,5); y = paramOP.binRMSE(:,2); ind = ~isnan(x);
    mdl2 = fitlm(x(ind), y(ind));
    gridOP.mdl2PVal(param1Ind, param2Ind) = mdl2.Coefficients.pValue("x1");
    gridOP.mdl2Slope(param1Ind, param2Ind) = mdl2.Coefficients.Estimate("x1");
    if (isequal(analysisParam.paramPoint, [param1Ind, param2Ind]) || isequal(analysisParam.paramPoint, "all")) ...
            && analysisParam.plotExtraFigs
        switch analysisParam.figSettings
            case 'manuscript'
                myPlotSettings(width=3.0, height=2.5)
            otherwise
                myPlotSettings
        end
        % [xBin, binRMSE', isActiveOnsetBin, fracBinDtWithActiveCluster', fracBinDtWithActiveClusterExpanded', binNormRMSE']
        binRMSELabels = {'decodeBin', 'Bin RMSE', 'isActiveOnsetBin', 'Frac. bin with active cluster ', 'Frac. bin with active cluster (most-active)', 'Bin RMSE (norm)'};
        xInd = 4;
        yInd = 2;
        x = paramOP.binRMSE(:, xInd);
        y = paramOP.binRMSE(:, yInd);
        ind = ~isnan(x); % [x~=0 & x~=1];
        mdl = fitlm(x(ind), y(ind));
        disp(mdl);
        figure;
        % plot(mdl); legend off
        scatter(x, y, 7.5)
        hline = refline(mdl.Coefficients.Estimate(2),  mdl.Coefficients.Estimate(1));
        hline.Color = 'k';
        hline.LineStyle = "--";
        xlabel(binRMSELabels{xInd})
        ylabel(binRMSELabels{yInd})
        title(...
            { ...
            ['intercept p=', num2str(mdl.Coefficients.pValue(1))], ...
            ['slope p=', num2str(mdl.Coefficients.pValue(2))] ...
            }, ...
            'FontWeight', 'normal', ...
            'FontSize', 12 ...
            )
        keyboard

        %{
        mdl = fitlm(adjacentBinFracClustActivePreceding, adjacentBinSimilaritiesCircShiftNorm);
        disp(mdl)
        figure;
        plot(mdl)
        legend off
        xlabel('Active cluster frac. (preceding bin)')
        ylabel('Adj. bin EMD, circ. adjusted (norm.)')
        %}

    end

    % Plot and analyze spikes/cell/PBE data in op struct
    grid.updateEventSpikingPlots( ...
        paramGridFigHandles, analysisParam, gridOP, paramOP, simParam ...
        )
end

runTime = toc;
disp(['Runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS')])

%% Plots analyzing hypothesis for prepaly

if analysisParam.plotExtraFigs
    myPlotSettings(width=3, height=2)
    xVals = simParam.variedParam(1).range;
    xLabel = simParam.variedParam(1).name;
    yVals = simParam.variedParam(2).range;
    yLabel = simParam.variedParam(2).name;
    nanMask = isnan(gridOP.nClustActive);

    figure;
    imagesc(xVals, yVals, gridOP.mdl1PVal', 'alphadata', ~nanMask')
    title('mdl1PVal')
    clim([0, 0.1]);
    colorbar;
    xlabel(xLabel); ylabel(yLabel)
    set(gca,'YDir','normal');

    figure;
    imagesc(xVals, yVals, gridOP.mdl1Slope', 'alphadata', ~nanMask')
    title('mdl1Slope')
    colorbar;
    xlabel(xLabel); ylabel(yLabel)
    set(gca,'YDir','normal');

    figure;
    imagesc(xVals, yVals, gridOP.mdl2PVal', 'alphadata', ~nanMask')
    title('mdl2PVal')
    clim([0, 0.1]);
    colorbar;
    xlabel(xLabel); ylabel(yLabel)
    set(gca,'YDir','normal');

    figure;
    imagesc(xVals, yVals, gridOP.mdl2Slope', 'alphadata', ~nanMask')
    title('mdl2Slope')
    colorbar;
    xlabel(xLabel); ylabel(yLabel)
    set(gca,'YDir','normal');

    figure;
    imagesc(xVals, yVals, gridOP.coactiveSymmetryPVal', 'alphadata', ~nanMask')
    title('coactiveSymmetryPVal')
    clim([0, 0.1]);
    colorbar;
    xlabel(xLabel); ylabel(yLabel)
    set(gca,'YDir','normal');

    figure;
    imagesc(xVals, yVals, gridOP.coactiveSymmetryIndex', 'alphadata', ~nanMask')
    title('coactiveSymmetryIndex')
    colorbar;
    xlabel(xLabel); ylabel(yLabel)
    set(gca,'YDir','normal');
end
