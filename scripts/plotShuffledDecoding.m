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

addpath(fullfile("..", "src"))
Config = utils.setConfig;
dataPath = Config.simDataPath;

filenameSplit = split(variedDecodeFilename, 'shuffledDecode');
gridFilename = filenameSplit{1};

load(fullfile(dataPath, gridFilename, variedDecodeFilename))

% Original, reference decoding result
tmp = dir(fullfile(dataPath, gridFilename, [gridFilename, 'decode_*']));
referenceDecodeFilename = fullfile(tmp.folder, tmp.name);
refDecode = load(referenceDecodeFilename);


%% Analysis options

% Analysis
analysisParam.ithEnv = 1;
analysisParam.nCIBootstraps = 1000;
analysisParam.bootFun = @median;

% Plotting
analysisParam.figSettings = 'manuscript';
analysisParam.plotEachCDF = false;
analysisParam.plotOldFigs = false;
analysisParam.plotExtraFigs = false;

%% Set up

switch analysisParam.figSettings
    case 'manuscript'
        myPlotSettings(width=1.75, height=1.5)
    otherwise
        myPlotSettings(width=2.5, height=2)
end

% Initialize output struct
op.pval = nan(1, analysisParam.nReplicates);
op.medianShift = nan(1, analysisParam.nReplicates);
op.medianShiftNorm = nan(1, analysisParam.nReplicates);
op.eventCount = nan(1, analysisParam.nReplicates);

disp(['Shuffle method: ', analysisParam.shuffleMethod])

%% Loop over decode parameters and collect statistics

for ithReplicate = 1:analysisParam.nReplicates

    allNet_rs = [];
    allNet_rs_shuffle = [];
    for ithNet = 1:modelParam.nNets
        if isempty(resultsStruct_variedDecode(ithReplicate).net{ithNet}{1}.replaytrajectory.weightedR2)
            continue % Skip this net if no events
        end
        net_rs = sqrt(abs(resultsStruct_variedDecode(ithReplicate).net{ithNet}{1}.replaytrajectory.weightedR2(:,analysisParam.ithEnv)));
        tmp1 = vertcat(resultsStruct_variedDecode(ithReplicate).net{ithNet}{1}.replaytrajectory.shuffle_weightedR2{:});
        tmp2 = tmp1(:,analysisParam.ithEnv);
        net_rs_shuffle = sqrt(abs([tmp2{:}]));
        allNet_rs = [allNet_rs, net_rs'];
        allNet_rs_shuffle = [allNet_rs_shuffle, net_rs_shuffle];
        if 0
            figure; hold on; ecdf(net_rs); ecdf(net_rs_shuffle)
        end
    end

    medianShift = median(allNet_rs)-median(allNet_rs_shuffle);
    medianShiftNorm = medianShift/median(allNet_rs_shuffle);
    [~,P,~] = kstest2(allNet_rs, allNet_rs_shuffle);

    op.eventCount(ithReplicate) =  numel(allNet_rs);
    op.pval(ithReplicate) = P;
    op.medianShift(ithReplicate) = medianShift;
    op.medianShiftNorm(ithReplicate) = medianShiftNorm;

    if analysisParam.plotEachCDF
        figure; hold on; title(['ithRep=', num2str(ithReplicate), ' pval=', num2str(P)])
        ecdf(allNet_rs); ecdf(allNet_rs_shuffle)
    end
end


%% Reference decoding result

%ithNet = 1;
%refDecode.resultsStruct(analysisParam.ithSimParam1, analysisParam.ithSimParam2, ithNet).results.replaytrajectory

allNet_rs = [];
allNet_rs_shuffle = [];
for ithNet = 1:modelParam.nNets
    if isempty(refDecode.resultsStruct(analysisParam.ithSimParam1, analysisParam.ithSimParam2, ithNet).results.replaytrajectory.weightedR2)
        continue % Skip this net if no events
    end
    net_rs = sqrt(abs(refDecode.resultsStruct(analysisParam.ithSimParam1, analysisParam.ithSimParam2, ithNet).results.replaytrajectory.weightedR2(:,analysisParam.ithEnv)));
    tmp1 = vertcat(refDecode.resultsStruct(analysisParam.ithSimParam1, analysisParam.ithSimParam2, ithNet).results.replaytrajectory.shuffle_weightedR2{:});
    tmp2 = tmp1(:,analysisParam.ithEnv);
    net_rs_shuffle = sqrt(abs([tmp2{:}]));
    allNet_rs = [allNet_rs, net_rs'];
    allNet_rs_shuffle = [allNet_rs_shuffle, net_rs_shuffle];
    if 0
        figure; hold on; ecdf(net_rs); ecdf(net_rs_shuffle)
    end
end

reference.medianShift = median(allNet_rs)-median(allNet_rs_shuffle);
reference.medianShiftNorm = reference.medianShift/median(allNet_rs_shuffle);
[~,reference.pval,KSSTAT] = kstest2(allNet_rs, allNet_rs_shuffle);
reference.eventCount =  numel(allNet_rs);


%% Plot

% Set up
yLabelString = 'Shuffles (count)';


% P-value of KS-test
y = op.pval';
refVal = reference.pval;
xLabelString = "Preplay p-value";
ci = bootci(analysisParam.nCIBootstraps, analysisParam.bootFun, y);
shuffPrctile = mean(reference.pval>y);
figure;
hold on
if max(y)>0.2
    histBinning = 0:0.25:1.0;
else
    histBinning = ceil(sqrt(numel(y)));
    y = log10(y);
    refVal = log10(refVal);
    xLabelString = "Preplay p-value (log10)";
end
h = histogram(y, histBinning);
xline(refVal, 'r:', 'lineWidth', 2)
ylabel(yLabelString)
xlabel(xLabelString);
title({ ...
    ['Actual: ', num2str(reference.pval)],...
    ['Shuffle range ', num2str(min(y)), ' to ', num2str(max(y))] ...
    ['Percentile of shuffle ', num2str(shuffPrctile)]
    }, ...
    'FontWeight', 'normal', ...
    'FontSize', 11 ...
    )

% Shift in median decode correlation
y = op.medianShift';
ci = bootci(analysisParam.nCIBootstraps, analysisParam.bootFun, y);
shuffPrctile = mean(reference.medianShift>y);
figure;
hold on
histogram(y);
xline(reference.medianShift, 'r:', 'lineWidth', 2)
ylabel(yLabelString)
xlabel("Preplay median r shift");
title({ ...
    ['Actual: ', num2str(reference.medianShift, 4)],...
    ['Shuffle range ', num2str(min(y)), ' to ', num2str(max(y))] ...
    ['Percentile of shuffle ', num2str(shuffPrctile)]
    }, ...
    'FontWeight', 'normal', ...
    'FontSize', 11 ...
    )


%% Calculate statistics
if isequal(shuffleToPlot, "clusterIndependent")

    % Check if p-value distribution is consistent with uniform dist [0, 1]
    x = op.pval';
    %Create expected counts variable, which is repeated 4 times as that is the
    %number of discrete values
    expCounts = repmat(numel(x)/4,[4,1]);
    %run the chi2gof test
    [h, chi2_p, chi2_stats]=chi2gof(x,'Ctrs',[0.125:0.25:1.0],'Expected',expCounts)
    disp(['p-val dist chi2gof, p=', num2str(chi2_p), ', chi2stat=', num2str(chi2_stats.chi2stat)])

    % Check if median shift mean is consistent with 0
    [h, p, shift_ci, stats] = ttest(op.medianShift);
    disp(['t-test 95% CI of median shift: ', num2str(shift_ci)])

elseif isequal(shuffleToPlot, "withinCluster")

    disp(['Out of ', num2str(analysisParam.nReplicates), ' shuffles, '])
    fracPValSig = mean(op.pval<0.05);
    disp(['Fraction of shuffle p-values that are significant: ', num2str(fracPValSig)])

    actualPerctOfShuff = mean(reference.medianShift>op.medianShift);
    disp(['Fraction of shuffles'' median shift less than actual: ', num2str(actualPerctOfShuff)])


elseif isequal(shuffleToPlot, "singleClusterCells")

    % Check if all shuffle p-values are significant
    disp(['Out of ', num2str(analysisParam.nReplicates), ' shuffles, '])
    fracPValSig = mean(op.pval<0.05);
    disp(['Fraction of shuffle p-values that are significant: ', num2str(fracPValSig)])
    disp('')

    % Check if median shift mean is consistent with reference value
    [h, p, shift_ci, stats] = ttest(op.medianShift-reference.medianShift);
    disp('t-test of median shift relative to reference value:')
    disp(['p=', num2str(p), ', 95% CI: ', num2str(shift_ci)])

end


%% Plot extra figures
if analysisParam.plotExtraFigs

    % Number of events included in decoding
    y = op.eventCount';
    ci = bootci(analysisParam.nCIBootstraps, analysisParam.bootFun, y);
    figure;
    hold on
    histogram(y);
    xline(reference.eventCount, 'r:', 'lineWidth', 2)
    ylabel(yLabelString)
    xlabel("Number of events");
    title({ ...
        ['Actual: ', num2str(reference.eventCount)],...
        ['Shuffle CI: ', num2str(ci(1)), ', ', num2str(ci(2)) ], ...
        ['Min: ', num2str(min(y)), ', Max: ', num2str(max(y))] ...
        }, ...
        'FontWeight', 'normal', ...
        'FontSize', 11 ...
        )

    % Normalized shift in median decode correlation
    y = op.medianShiftNorm';
    ci = bootci(analysisParam.nCIBootstraps, analysisParam.bootFun, y);
    shuffPrctile = mean(reference.medianShiftNorm>y);
    figure;
    hold on
    histogram(y);
    xline(reference.medianShiftNorm, 'r:', 'lineWidth', 2)
    ylabel(yLabelString)
    xlabel("Preplay median r shift (norm.)");
    title({ ...
        ['Actual (red): ', num2str(reference.medianShiftNorm, 4)],...
        ['Shuffle range ', num2str(min(y)), ' to ', num2str(max(y))] ...
        ['Percentile of shuffle ', num2str(shuffPrctile)]
        }, ...
        'FontWeight', 'normal', ...
        'FontSize', 11 ...
        )
end

