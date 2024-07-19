%% Plot results of varyDecodeParam.m
% plotVariedDecodeParam.m
%


%% Choose which varied decode parameter result to analyze

% Creates Figure 4—figure supplement 2
% Decoding "mainPaperGrid" with 25 random subset replicates at each point
variedDecodeFilename = 'grid_2023-07-22T14-50variedDecode2023-09-10T21-04.mat';


%% Load results

addpath(fullfile("..", "src"))
Config = utils.setConfig;
dataPath = Config.simDataPath;

filenameSplit = split(variedDecodeFilename, 'variedDecode');
gridFileName = filenameSplit{1};

load(fullfile(dataPath, gridFileName, variedDecodeFilename))


%% Set up

analysisParam.useCellCountLabel = true;
analysisParam.figSettings = 'manuscript';
analysisParam.plotEachCDF = false;
analysisParam.plotOldFigs = false;
analysisParam.ithEnv = 1;
analysisParam.nCIBootstraps = 1000;
analysisParam.bootFun = @median;
analysisParam.markerSize = 20;

switch analysisParam.figSettings
    case 'manuscript'
        myPlotSettings(width=2, height=1.75)
    otherwise
        myPlotSettings(width=2.5, height=2)
end

op.pval = nan(numel(analysisParam.variedParam.range), analysisParam.nReplicates);
op.medianShift = nan(numel(analysisParam.variedParam.range), analysisParam.nReplicates);
op.medianShiftNorm = nan(numel(analysisParam.variedParam.range), analysisParam.nReplicates);
op.eventCount = nan(numel(analysisParam.variedParam.range), analysisParam.nReplicates);


%% Loop over decode parameters and collect statistics

for ithParam = 1:numel(analysisParam.variedParam.range)
    for ithReplicate = 1:analysisParam.nReplicates

        allNet_rs = [];
        allNet_rs_shuffle = [];
        for ithNet = 1:modelParam.nNets
            if isempty(resultsStruct_variedDecode(ithParam,ithReplicate).net{ithNet}{1}.replaytrajectory.weightedR2)
                continue % Skip this net if no events
            end
            net_rs = sqrt(abs(resultsStruct_variedDecode(ithParam,ithReplicate).net{ithNet}{1}.replaytrajectory.weightedR2(:,analysisParam.ithEnv)));
            tmp1 = vertcat(resultsStruct_variedDecode(ithParam,ithReplicate).net{ithNet}{1}.replaytrajectory.shuffle_weightedR2{:});
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
        [~,P,KSSTAT] = kstest2(allNet_rs, allNet_rs_shuffle);

        op.eventCount(ithParam, ithReplicate) =  numel(allNet_rs);
        op.pval(ithParam, ithReplicate) = P;
        op.medianShift(ithParam, ithReplicate) = medianShift;
        op.medianShiftNorm(ithParam, ithReplicate) = medianShiftNorm;

        if analysisParam.plotEachCDF
            figure; hold on; title(['ithParam=', num2str(ithParam), ', ithRep=', num2str(ithReplicate), ' pval=', num2str(P)])
            ecdf(allNet_rs); ecdf(allNet_rs_shuffle)
        end
    end
end


%% Plot

% Setup for plotting
xVals = analysisParam.variedParam.range;
if analysisParam.useCellCountLabel
    xlabelString = 'Cells (count)';
    xVals = xVals.*modelParam.n_E;
else
    xlabelString = 'Fraction of cells';
end
X = repmat(xVals, analysisParam.nReplicates, 1);
%xtickvals = [[12, 25, 50, 100, 200, 375]];
xtickvals = [sort(round(xVals([2:3:end]) ))]; % sort(round(xVals))
%xtickvals = [sort(round(xVals([end:-5:1]) ))]; % sort(round(xVals))

% Number of events included in decoding
y = op.eventCount';
ci = bootci(analysisParam.nCIBootstraps,analysisParam.bootFun,y);
figure; hold on
scatter(X, y, analysisParam.markerSize, 'r');
%errorbar(xVals, analysisParam.bootFun(y), ci(1,:)-analysisParam.bootFun(y), ci(2,:)-analysisParam.bootFun(y), 'k');
plot(xVals, analysisParam.bootFun(y), 'k')
xlabel(xlabelString); ylabel('Number of events')
%set(gca, 'xdir', 'reverse')
set(gca,'XScale', 'Log')
xticks(xtickvals); xtickangle(45)

% P-value of KS-test
y = op.pval';
ci = bootci(analysisParam.nCIBootstraps,analysisParam.bootFun,y);
figure; hold on
scatter(X, y, analysisParam.markerSize, 'r');
%errorbar(xVals, analysisParam.bootFun(y), ci(1,:)-analysisParam.bootFun(y), ci(2,:)-analysisParam.bootFun(y), 'k');
plot(xVals, analysisParam.bootFun(y), 'k')
xlabel(xlabelString); ylabel('Preplay p-value')
yline(0.05, ':', '0.05', LabelHorizontalAlignment='right')
set(gca,'YScale', 'Log')
set(gca,'XScale', 'Log')
xticks(xtickvals); xtickangle(45)

% Shift in median decode correlation
y = op.medianShift';
ci = bootci(analysisParam.nCIBootstraps,analysisParam.bootFun,y);
figure; hold on
scatter(X, y, analysisParam.markerSize, 'r');
%errorbar(xVals, analysisParam.bootFun(y), ci(1,:)-analysisParam.bootFun(y), ci(2,:)-analysisParam.bootFun(y), 'k');
plot(xVals, analysisParam.bootFun(y), 'k')
xlabel(xlabelString); ylabel('Preplay median r shift')
yline(0.0, ':')
set(gca,'XScale', 'Log')
xticks(xtickvals); xtickangle(45)

% Number of subsets with p<0.05
y = op.pval';
figure; plot(X, sum(y<0.05), '-x')
xlabel(xlabelString); ylabel('subsets with p<0.05 (count)')
set(gca,'XScale', 'Log')
xticks(xtickvals); xtickangle(45)


%%  Old figure versions, not used in mansucript

if analysisParam.plotOldFigs
    y = op.eventCount';
    figure; hold on
    scatter(X, y, analysisParam.markerSize, 'r');
    %errorbar(xVals, mean(y), std(y), 'k')
    plot(xVals, mean(y), 'k')
    xlabel(xlabelString); ylabel('Number of events')
    %set(gca, 'xdir', 'reverse')

    y = op.pval';
    figure; hold on
    scatter(X, y, analysisParam.markerSize, 'r');
    %h = errorbar(xVals, mean(y), std(y), 'k');
    %ylim manual
    %h.LData = h.YData - max(eps/10,h.YData-h.LData); % for log of negative data
    plot(xVals, mean(y), 'k')
    xlabel(xlabelString); ylabel('Preplay p-value')
    yline(0.05, ':', '0.05', LabelHorizontalAlignment='right')
    set(gca,'YScale', 'Log')
    %set(gca, 'xdir', 'reverse')

    y = op.medianShift';
    figure; hold on
    scatter(X, y, analysisParam.markerSize, 'r');
    %errorbar(xVals, mean(y), std(y), 'k')
    plot(xVals, mean(y), 'k')
    xlabel(xlabelString); ylabel('Preplay median r shift')
    yline(0.0, ':')
    %set(gca, 'xdir', 'reverse')

    y = op.medianShiftNorm';
    figure; hold on
    scatter(X, y, analysisParam.markerSize, 'r');
    errorbar(xVals, mean(y), std(y), 'k')
    xlabel(xlabelString); ylabel('Decode r median shift (norm.)')
    yline(0.0, ':')
    %set(gca, 'xdir', 'reverse')
end
