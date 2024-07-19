% Plots preplay decoding results from shuffling cell identities
%
% The script runShuffledDecoding.m generates the files this script plots.


%% Choose which varied decode parameter result to analyze

shuffleToPlot = "singleClusterCells";


%% Load results

switch shuffleToPlot
    case "clusterIndependent"
        variedDecodeFilename = 'grid_2023-07-22T14-50shuffledDecode2024-04-26T21-15';
    case "withinCluster"
        variedDecodeFilename = 'grid_2023-07-22T14-50shuffledDecode2024-04-29T04-21';
    case "singleClusterCells"
        variedDecodeFilename = 'grid_2023-07-22T14-50shuffledDecode2024-05-01T05-18';
    otherwise
        error("Unknown value of shuffleToPlot")
end


%% Load results

addpath(fullfile("..", "src"))
Config = utils.setConfig;
dataPath = Config.simDataPath;

filenameSplit = split(variedDecodeFilename, 'shuffledDecode');
gridFilename = filenameSplit{1};

%
load(fullfile(dataPath, gridFilename, variedDecodeFilename))

% Original, reference decoding result
tmp = dir(fullfile(dataPath, gridFilename, [gridFilename, 'decode_*']));
%referenceDecodeFilename = fullfile(tmp.folder, tmp.name);
%refDecode = load(referenceDecodeFilename);

% Load reference decoding parameters and results
if exist('shuffledDecodingReferenceDecodeName', 'var')
    disp('Assuming loaded reference decode is correct')
else
    tmpDir = dir(fullfile(dataPath, gridFilename, [gridFilename, 'decode_*']));
    shuffledDecodingReferenceDecodeName = tmpDir.name(1:end-4);
    [refDecode.modelParam, ...
        refDecode.simParam, ...
        refDecode.resultsStruct, ...
        refDecode.PFresultsStruct, ...
        refDecode.preplayGridSpikes, ...
        ] = grid.loadResults(shuffledDecodingReferenceDecodeName);
end


%% Analysis options

% Analysis
analysisParam.ithEnv = 1;
analysisParam.ithPreplayTrial = 1;
analysisParam.nCIBootstraps = 1000;
analysisParam.bootFun = @median;

% Event spike rank analyses
analysisParam.minPfRate = 3; % Exclude cells with PF peaks below this rate
analysisParam.pfRankMethod = 'mean'; % peak or mean
analysisParam.invertDirectionality = true; % if true, invert time in events with negative decode slopes
analysisParam.onlySingleClustCells = false; % if true, exclude cells that belong to multiple clusters from the analysis

analysisParam.matchPFtoShuffle = false; % if true, shuffle PFs to match the spike shuffling
analysisParam.usePfCmUnits = true; % Use cm if true, otherwise spatial bin

analysisParam.runRasterChecks = false;
analysisParam.calcCorrForEachNet = false;
analysisParam.plotExtraFigs = false;

% Plotting
analysisParam.figSettings = 'manuscript';
analysisParam.plotEachCDF = false;
analysisParam.plotOldFigs = false;


%% Set up

switch analysisParam.figSettings
    case 'manuscript'
        myPlotSettings(width=2.0, height=1.75)
    otherwise
        myPlotSettings(width=2.5, height=2)
end

% Set up indices
param1Ind = analysisParam.ithSimParam1;
param1Val = refDecode.simParam.variedParam(1).range(param1Ind);
param1Name = refDecode.simParam.variedParam(1).name;
param2Ind = analysisParam.ithSimParam2;
param2Val = refDecode.simParam.variedParam(2).range(param2Ind);
param2Name = refDecode.simParam.variedParam(2).name;
linearParamInd = find( ...
    all([param1Val; param2Val] == refDecode.simParam.parameterSets_vec, 1) ...
    );

% Set up modelParam
paramPointModelParam = modelParam;
paramPointModelParam.(param1Name) = param1Val;
paramPointModelParam.(param2Name) = param2Val;
paramPointModelParam = set_depedent_parameters(paramPointModelParam);

% Initialize output struct
op = struct;

disp(['Shuffle method: ', analysisParam.shuffleMethod])


%% Loop over decode parameters and collect statistics

for ithReplicate = 1:analysisParam.nReplicates

    % Match spike shuffling
    paramPointSpikes = refDecode.preplayGridSpikes{linearParamInd};

    % Need to re-format to expected form and shuffle cell ids in spikes
    tmpResults = struct;
    tmpPFResults = struct;
    replicateShuffledSpikes = paramPointSpikes;
    for ithNet = 1:modelParam.nNets
        netPreplayResults = resultsStruct_variedDecode(ithReplicate).net{ithNet}{1};
        netPfResults = PFresultsStruct_variedDecode(ithReplicate).net{ithNet}{1};
        netSpikes = squeeze(paramPointSpikes(ithNet,analysisParam.ithPreplayTrial,:));

        % Shuffle spikes using saved ordering from decode_events.m
        cellidxm = netPreplayResults.replaytrajectory.cellidxm;
        shuffledCellidxm = netPreplayResults.replaytrajectory.shuffledCellidxm;
        hpnum = size(shuffledCellidxm, 1);
        day = 1;
        epoch = 1;
        shuffledNetSpikes = netSpikes;
        for ithCell = 1:hpnum
            shuffledNetSpikes{cellidxm(ithCell,2)} = ...
                netSpikes{shuffledCellidxm(ithCell,2)};
        end

        % Save results in structure needed for +grid functions
        tmpResults(ithNet).results = netPreplayResults;
        tmpPFResults(ithNet).results{1} = netPfResults;
        replicateShuffledSpikes(ithNet, analysisParam.ithPreplayTrial,:) = shuffledNetSpikes;
    end

    if ithReplicate == 1
        plotFigs=true;
    else
        plotFigs=false;
    end
    replicateOutput = grid.calculateEventSpikingSequences(...
        analysisParam, ...
        paramPointModelParam, ...
        tmpResults, ...
        tmpPFResults, ...
        replicateShuffledSpikes, ...
        plotFigs=plotFigs, ...
        matchPFtoShuffle=analysisParam.matchPFtoShuffle ...
        );

    op.netMdlPval(ithReplicate) = replicateOutput.netMdl.Coefficients.pValue(2);
    op.clustMdlPval(ithReplicate) = replicateOutput.clustMdl.Coefficients.pValue(2);
    op.binMdlPval(ithReplicate) = replicateOutput.binMdl.Coefficients.pValue(2);

    op.netMdlSlope(ithReplicate) = replicateOutput.netMdl.Coefficients.Estimate(2);
    op.clustMdlSlope(ithReplicate) = replicateOutput.clustMdl.Coefficients.Estimate(2);
    op.binMdlSlope(ithReplicate) = replicateOutput.binMdl.Coefficients.Estimate(2);

    disp(['Replicate ', num2str(ithReplicate), ' completed'])
end


%% Reference decoding result

% Calculate event spike sequence stats
referenceOutput = grid.calculateEventSpikingSequences(...
    analysisParam, ...
    paramPointModelParam, ...
    refDecode.resultsStruct(param1Ind, param2Ind, :), ...
    refDecode.PFresultsStruct(param1Ind, param2Ind, :), ...
    refDecode.preplayGridSpikes{linearParamInd} ...
    );


%% Plot


% Set up
yLabelString = 'Shuffles (count)';
xLabelString = 'Regression slope';

y = op.netMdlSlope;
refVal = referenceOutput.netMdl.Coefficients.Estimate(2);
plotStatHist(y, refVal)
ylabel(yLabelString);
%xlabel(xLabelString);
xlabel({'Within-network', xLabelString})
%figure; hold on; ecdf(y); xline(refVal, 'r:', 'lineWidth', 2)
[~,p] = ttest(y);
disp(['netMdl: refVal=', num2str(refVal), ', p=', num2str(p)])

y = op.clustMdlSlope;
refVal = referenceOutput.clustMdl.Coefficients.Estimate(2);
plotStatHist(y, refVal)
ylabel(yLabelString);
%xlabel(xLabelString);
xlabel({'Within-cluster', xLabelString})
%figure; hold on; ecdf(y); xline(refVal, 'r:', 'lineWidth', 2)
[~,p] = ttest(y);
disp(['clustMdl: refVal=', num2str(refVal), ', p=', num2str(p)])

y = op.binMdlSlope;
refVal = referenceOutput.binMdl.Coefficients.Estimate(2);
plotStatHist(y, refVal)
ylabel(yLabelString);
%xlabel(xLabelString);
xlabel({'Within-bin', xLabelString})
%figure; hold on; ecdf(y); xline(refVal, 'r:', 'lineWidth', 2)
[~,p] = ttest(y);
disp(['binMdl: refVal=', num2str(refVal), ', p=', num2str(p)])

function figHandle = plotStatHist(y, refVal)
shuffPrctile = mean(refVal>y);
figure;
hold on
histogram(y);
%histogram(y, 0:0.25:1.0);
xline(refVal, 'r:', 'lineWidth', 2)
%ylabel(yLabelString)
xlabel("Number of events");
title({ ...
    ['Actual: ', num2str(refVal)],...
    ['Shuffle range ', num2str(min(y)), ' to ', num2str(max(y))] ...
    ['Percentile of shuffle ', num2str(shuffPrctile)]
    }, ...
    'FontWeight', 'normal', ...
    'FontSize', 11 ...
    )
%xticklabels(xticks)
end