%% Plot preplay decoding statistics across grid search
% gridAnalysisDecoding.m
%

%% Choose which simulation results to analyze

decodeFileID = "mainPaperGrid"; % Select which decoding file result ID to analyze

figSettings = 'manuscript0'; % manuscript0 for Figure 4/6, manuscript1 for figure R5


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
    PFresultsStruct ...
    ] = grid.loadResults(decodeName);


%% Analysis parameters
analysisParam.paramPoint = [2, 3]; % Select which example param point to plot
analysisParam.ithEnv        = 1; % which environments decoding and place fields to plot/analyze

% Plotting options
analysisParam.figSettings = figSettings;
analysisParam.MSPlots = true;
analysisParam.plotExtraFigs = true; % Only plot manuscript figures if false
analysisParam.plotAllCDFs   = false;
analysisParam.plotCombinedCDFs = false;

% Analysis options
analysisParam.maxNEvents    = inf; % downsample number of replay events for kstest to maxNEvents, use inf for no downsampling
analysisParam.combineNetData = true;
analysisParam.useWeightedDecode = true; % slope and R2 for correlations by either peak prob or weighted prob
analysisParam.SWIcorrYClass     = 'KSmedian'; % 'pval', 'KSmedian', or 'fracSig'
analysisParam.corrLogVals       = false; % Can't be true with SWIcorrYClass='KSmedian'
analysisParam.minNEventsSWI = 10;
analysisParam.excludeFullConn   = true;
analysisParam.fullConnBoundary  = @(mnc) (mnc.^2)/modelParam.conn_prob;
analysisParam.compPairWiseEnv = true; % Option for when comparing decodes from multiple environments
analysisParam.capacitySignificanceAlpha = 0.01; %0.05/2;

disp(analysisParam)

assert(any(strcmp(analysisParam.figSettings, {'standard', 'manuscript0', 'manuscript1', 'SfNPoster', 'manuscriptAlt'})))


%% Initialize the output structure op_comb with desired fields
%   op_comb.KSpval(ithParam1, ithParam2) % combined data at this param point
PFdecode = {'KSpval', 'KSstat', 'fracEventSig', 'KSMedianDiff', 'multiEnvPval',  ...
    'envComp_p', 'envComp_shift', 'envCapacityPred'};
opCombFields = [PFdecode, {'other'}];

tempCell = cell(length(opCombFields),1);
tempCell(:)= {nan(numel(simParam.variedParam(1).range), numel(simParam.variedParam(2).range))};
opComb = cell2struct(tempCell,opCombFields);


%% Initialize the output structure op with desired fields
%   op.SWI(ithParam1, ithParam2, ithNet)
decodeTest = {'KSpval', 'KSstat', 'KSMedianDiff', 'netGroupID'};
decodeStats = {'decodeSlope', 'decodeSlopeZScore', 'decodeEntropy', 'decodeVar'};
opFields = [decodeTest, decodeStats, {'SWI', 'fracEventSig', 'nEvents'}];

tempCell = cell(length(opFields),1);
tempCell(:)= {nan(numel(simParam.variedParam(1).range), numel(simParam.variedParam(2).range), simParam.nNets)};
op = cell2struct(tempCell,opFields);


%% Loop over all parameter sets and all networks

% Main analysis function
rng('default')
tic
[op, opComb] = loopoOPCalc(op, opComb, analysisParam, simParam, modelParam, PFresultsStruct, resultsStruct);
runTime = toc;
disp(['Loop runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS')])

% TEMP, to fix in original simulateParameterGrid creation of it
% Fill missing values in resultsStruct with nan
fieldNames = fieldnames(resultsStruct(1,1,1).results);

for ithLinInd = 1:numel(resultsStruct)
    for ithFieldName = 1:length(fieldNames)
        thisFieldName = fieldNames{ithFieldName};
        if ~isfield(resultsStruct(ithLinInd).results, thisFieldName)
            resultsStruct(ithLinInd).results.(thisFieldName)=nan;
        end
    end
end

% Some additional preplay statistics
preplayStats.frac_participation = arrayfun(@(x) x.results.frac_participation, resultsStruct); % mean fraction of cells firing per event
preplayStats.meanRate = arrayfun(@(x) x.results.meanRate, resultsStruct); % mean E-cell firing rate
preplayStats.eventLength = arrayfun(@(x) mean(x.results.eventLength), resultsStruct); % mean event length
preplayStats.eventFreq = arrayfun(@(x) x.results.numEvents/modelParam.t_max_preplay, resultsStruct); % event frequency of occurance

% Calculate some means over nets
opComb.decodeSlopeZScore = nanmean(op.decodeSlopeZScore, 3);
opComb.decodeEntropy = nanmean(op.decodeEntropy, 3);


%% Set up for plotting

xParamvec = simParam.variedParam(1).range;
xName = simParam.variedParam(1).name;
yParamvec = simParam.variedParam(2).range;
yName = simParam.variedParam(2).name;

if contains( xName , 'del_G_syn_' )
    xName = ['W_', xName(11:end), ' (pS)'];
    xParamvec = xParamvec*10^12;
end
if contains( yName , 'del_G_syn_' )
    yName = ['W_', yName(11:end), ' (pS)'];
    yParamvec = yParamvec*10^12;
end

if analysisParam.useWeightedDecode
    analysisTitle = 'KS-test, weighted decode corrs';
else
    analysisTitle = 'KS-test, peak pos. decode corrs';
end
if analysisParam.combineNetData
    cbLabel1 = 'Combined p-val';
    cbLabel2 = 'Combined KS-stat';
else
    cbLabel1 = 'Median p-val';
    cbLabel2 = 'Median KS-stat';
end

% Matrix of values to exclude if analysisParam.excludeFullConn
if analysisParam.excludeFullConn && isequal(simParam.variedParam(1).name, 'mnc') && isequal(simParam.variedParam(2).name, 'clusters')
    mncVec = simParam.variedParam(1).range;
    nClustVec = simParam.variedParam(2).range;
    excludedBadConnProb = [(nClustVec' > analysisParam.fullConnBoundary(mncVec))]';
else
    excludedBadConnProb = zeros(numel(simParam.variedParam(1).range), numel(simParam.variedParam(2).range));
end


%% Track capacity
if isfield(opComb, "envCapacityPred") && (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)

    if analysisParam.MSPlots
        myPlotSettings(width=2.25, height=1.5)
    end

    figure;
    imagesc(xParamvec, yParamvec, opComb.envCapacityPred', 'AlphaData', ~isnan(opComb.envCapacityPred')&~excludedBadConnProb')
    colorbar;
    set(gca,'YDir','normal');
    xlabel(xName,'Interpreter','none')
    ylabel(yName,'Interpreter','none')
    title("Predicted track capacity", 'FontWeight','Normal')
    % TODO: What is the expected null result? Does their analysis account
    % for the expected 5% false positive rate?
    %{
    k = 1:20;
    alpha = 0.05;
    y = alpha.^k;
    yy(1) = y(1);
    for i = 2:max(k)
        yy(i) = k(i)*y(1)-y(i);
    end
    figure
    plot(k, alpha*k)
    %}

end


%% Plot manuscript figures

if analysisParam.MSPlots

    % Event statistics across parameter grid
    switch analysisParam.figSettings
        case 'manuscript0'
            myPlotSettings(width=3.33, height=2.5, ttlfsz=1.0)
        case 'manuscript1'
            myPlotSettings(width=2.5, height=1.75, ttlfsz=1.0)
        otherwise
            myPlotSettings(width=2.5, height=2)
    end

    plotFields = {'frac_participation', 'meanRate', 'eventLength', 'eventFreq'};
    plotInds = [1:4];
    plotTitles = {'Frac. participation', 'Mean Rate (Hz)', 'Event Length', 'Event Freq.'};
    figure
    for ithField = 1:numel(plotFields)
        currentField = plotFields{ithField};
        ax = subplot(2,2,plotInds(ithField));
        imagesc(xParamvec, yParamvec, mean([preplayStats.(currentField)], 3)', 'AlphaData', ~isnan(mean([preplayStats.(currentField)], 3)' )&~excludedBadConnProb')
        set(gca,'YDir','normal')
        cb = colorbar(); cb.Label.String = '';
        %xlabel(xName,'Interpreter','none')
        %ylabel(yName,'Interpreter','none')
        title(plotTitles{ithField}, 'FontWeight','Normal')
        if any(strcmp(currentField, {'linearPeakRMSD', 'meanPFRMSD'} ))
            colormap(ax, flipud(parula))
        end
    end

    % Decode statistics across parameter grid
    % Modifications to fit to very small supplement panel for non-uniform
    % axis labels
    if numel(unique(diff(xParamvec)))==1
        xValVec = xParamvec;
        yValVec = yParamvec;
    else
        xValVec = 1:numel(xParamvec);
        yValVec = 1:numel(yParamvec);
    end

    if isequal(figSettings, 'manuscript1') && isequal(xParamvec(2), 2)
        xValVec = xValVec(1:end-1);
        yValVec = yValVec(1:end-1);
    end

    plotFields = {'KSpval', 'KSMedianDiff', 'fracEventSig', 'decodeEntropy'};
    plotInds = [1:4];
    oldPvalStyle = false;
    plotTitles = {'Decode p-value', 'Median shift', 'Frac. Significant', 'Decode Entropy'};
    figure
    for ithField = 1:numel(plotFields)
        currentField = plotFields{ithField};
        ax = subplot(2,2,plotInds(ithField));
        C = mean([opComb.(currentField)], 3)';
        bcp = excludedBadConnProb;
        if isequal(figSettings, 'manuscript1') && isequal(xParamvec(2), 2)
            C(simParam.variedParam(1).range==2,:) = [];
            C(:,simParam.variedParam(2).range==2) = [];
            bcp(simParam.variedParam(1).range==2,:) = [];
            bcp(:,simParam.variedParam(2).range==2) = [];
        end
        panelImage = imagesc(xValVec, yValVec, C, 'AlphaData', ~isnan(C)&~bcp');
        set(gca,'YDir','normal')
        cb = colorbar(); cb.Label.String = '';
        %xlabel(xName,'Interpreter','none')
        %ylabel(yName,'Interpreter','none')
        title(plotTitles{ithField}, 'FontWeight','Normal')
        if any(strcmp(currentField, {'linearPeakRMSD', 'meanPFRMSD', 'KSpval', 'decodeEntropy'} ))
            colormap(ax, flipud(parula))
        end
        if any(strcmp(currentField, {'KSpval'} ))
            if oldPvalStyle || isequal(figSettings, 'manuscript1')
                clim([0, 0.1])
            else
                set(gca,'ColorScale','log')
                % clim([min(panelImage.CData(:)), max(panelImage.CData(:))])
                clim([min(panelImage.CData(:)), 0.1])
                if isequal(cb.Ticks, 1e-50)
                    clim([1e-20, 0.1])
                    cb.Ticks = [1e-20, 1e-10, 1e-1];
                end
            end
        end
        if numel(unique(diff(xParamvec)))~=1
            tmp1 = xParamvec;
            tmp2 = yParamvec;
            if isequal(figSettings, 'manuscript1') && isequal(xParamvec(2), 2)
                tmp1(xParamvec==2) = [];
                tmp2(yParamvec==2) = [];
            end
            xticks(xValVec); xticklabels(tmp1)
            yticks(yValVec); yticklabels(tmp2)
        end
    end

end


%% Plot decode correlation stats:

if (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
    % myPlotSettings(width=3, height=1.5) % for poster
    % myPlotSettings(width=4, height=3, lw=2, afs=14, alw=2) % ppt format
    myPlotSettings(width=5, height=4)

    figure;
    sgtitle(analysisTitle)

    subplot(2,2,1)
    imagesc(xParamvec, yParamvec, opComb.KSpval', 'AlphaData', ~isnan(opComb.KSpval')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = cbLabel1;
    xlabel(xName,'Interpreter','none')
    ylabel(yName,'Interpreter','none')
    %title(analysisTitle)


    subplot(2,2,2)
    imagesc(xParamvec, yParamvec, opComb.KSstat', 'AlphaData', ~isnan(opComb.KSstat')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = cbLabel2;
    xlabel(xName,'Interpreter','none')
    ylabel(yName,'Interpreter','none')
    %title(analysisTitle)

    subplot(2,2,3)
    imagesc(xParamvec, yParamvec, opComb.fracEventSig', 'AlphaData', ~isnan(opComb.fracEventSig')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = 'Frac Sig';
    xlabel(xName,'Interpreter','none')
    ylabel(yName,'Interpreter','none')
    %title(analysisTitle)
    % caxis([prctile(op, 2.5, 'all'), prctile(op, 97.5, 'all')])
    % set(gca,'ColorScale','log')
    % hold on; plot(simParam.variedParam(1).range, exp(simParam.variedParam(1).range/1.1)-1); plot(simParam.variedParam(1).range, exp((simParam.variedParam(1).range-1)*5)+15);

    % Plot with better colormap
    logPvalData = log10( opComb.KSpval');

    ax = subplot(2,2,4);
    imagesc(xParamvec, yParamvec, logPvalData, 'AlphaData', ~isnan(opComb.KSpval')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = cbLabel2;
    xlabel(xName,'Interpreter','none')
    ylabel(yName,'Interpreter','none')
    title(analysisTitle)

    N = 256; n = N/2;
    cm = NaN(N,3);
    cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
    cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)'];
    cm(:,3) = [linspace(0,1,n)';ones(N-n,1)];
    alpha = 0.05;
    bonferroniCorr = false;
    if bonferroniCorr; alpha = 0.05./numel(logPvalData); end
    set(ax,'clim',[log10(alpha)*2 0])
    set(ax,'colormap',cm)
    cb = colorbar('Direction','reverse','Ticks',[log10(alpha/10),log10(alpha),log10(alpha*10)],'TickLabels',[alpha/10,alpha,alpha*10]);
    cb.Label.String = 'KS test p-value';
    title(''); xlabel('Mean cluster membership'); ylabel('Number of clusters')
    % set(gca, 'XTick',xParamvec, 'XTickLabel',num2str(xParamvec')) % For grid starting at mnc=1.5

    % ithYval = 1:numel(yParamvec);
    %  ithYval = [2, 4,  8,  12];
    ithYval = [1];
    xx = xParamvec;
    yy = opComb.KSpval';
    figure; hold on;
    plot(xx, yy, '-o'); yline(0.05)
    xlabel('Mean cluster membership'); ylabel('KS test p-value');
    set(gca, 'YScale', 'log')
    legend({num2str( yParamvec(ithYval)')}, 'Location', 'Best')
end

%% Extra parameter grid plots

if (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
    figure;
    sgtitle(analysisTitle)

    subplot(2,2,1)
    imagesc(xParamvec, yParamvec, nanmean(op.decodeSlope, 3)', 'AlphaData', ~isnan(nanmean(op.decodeSlope, 3)')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = 'Mean abs. decode slope';
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

    subplot(2,2,2)
    imagesc(xParamvec, yParamvec, nanmean(op.decodeSlopeZScore, 3)', 'AlphaData', ~isnan(nanmean(op.decodeSlopeZScore, 3)')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = 'Shuffle scored decode slope';
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

    subplot(2,2,3)
    imagesc(xParamvec, yParamvec, nanmean(op.decodeEntropy, 3)', 'AlphaData', ~isnan(nanmean(op.decodeEntropy, 3)')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = 'Mean entropy';
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

    subplot(2,2,4)
    imagesc(xParamvec, yParamvec, nanmean(op.decodeVar, 3)', 'AlphaData', ~isnan(nanmean(op.decodeVar, 3)')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = 'Mean decode variance';
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

    myPlotSettings
end

%% Plot net-wise scatter, if single parameter point was run above

if analysisParam.plotExtraFigs || ~analysisParam.MSPlots
    figure; scatterhist(op.KSMedianDiff(:), op.KSpval(:), 'Kernel','on') % , 'Group', op.netGroupID(:)
    title('All networks from all parameter points')
    xlabel('median(actual R^2)-median(shuff R^2)')
    ylabel('KStest p-value')

    figure; scatterhist(op.KSstat(:), op.KSpval(:), 'Kernel','on')
    title('All networks from all parameter points')
    xlabel('KS-statistic')
    ylabel('KStest p-value')

    try
        [a, b] = find(opComb.KSpval'<0.05); ind1 = [a]; ind2 = [b];
        ind1 = [2]; ind2 = [4];

        X = squeeze(op.KSMedianDiff(ind1,ind2,:));
        Y = squeeze(op.KSpval(ind1,ind2,:));
        ID = squeeze(op.netGroupID(ind1,ind2,:));
        figure; scatterhist(X(:), Y(:), 'Group', ID(:), 'Kernel','on')
        % figure; scatter(X(:), Y(:), [], ID(:))
        title(['Nets from parameter point index (', num2str(ind1), ', ', num2str(ind2), ')'])
        xlabel('median(actual R^2)-median(shuff R^2)')
        ylabel('KStest p-value')
    end
end

%% A few misc plots

if (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
    figure
    imagesc(xParamvec, yParamvec, opComb.KSMedianDiff', 'AlphaData', ~isnan(opComb.KSMedianDiff')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); %cb.Label.String = 'KSMedianDiff';
    title('KSMedianDiff');
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

    figure
    imagesc(xParamvec, yParamvec, squeeze(median(op.SWI, 3))', 'AlphaData', ~isnan(squeeze(median(op.SWI, 3))')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); %cb.Label.String = 'Mean decode variance';
    title('Small-world index (median)');
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')
end

%% SWI correlation analysis
%
% op.SWI: SWI for each network (3rd dimension)
% op.KSpval: decode cdf p-value for each network
% opComb.KSpval: decode p-value for each parameter point (networks combined)

disp(['Using minNEventsSWI=', num2str(analysisParam.minNEventsSWI)])

if analysisParam.MSPlots
    myPlotSettings(width=1.5, height=1.25); scatterMS = 15;
else
    myPlotSettings(width=3.0, height=2.0); scatterMS = 35;
end

mncInd = find(strcmp({simParam.variedParam.name}, 'mnc'));
nClustInd = find(strcmp({simParam.variedParam.name}, 'clusters'));
if ~isempty(mncInd) && ~isempty(nClustInd)
    mncVec = simParam.variedParam(mncInd).range;
    nClustVec = simParam.variedParam(nClustInd).range;

    if analysisParam.excludeFullConn
        SWIParamsToInclude = [~(nClustVec' > analysisParam.fullConnBoundary(mncVec))]';
    else
        SWIParamsToInclude = true(size(opComb.KSpval));
    end
    SWINetsToInclude = repmat(SWIParamsToInclude, 1, 1, simParam.nNets);

    SWINetsToInclude = SWINetsToInclude & (op.nEvents>=analysisParam.minNEventsSWI);
    SWIParamsToInclude = SWIParamsToInclude & (squeeze(nansum(op.nEvents, 3))>=analysisParam.minNEventsSWI);

    % Net-wise correlation:
    X_allNets = op.SWI(SWINetsToInclude);
    switch analysisParam.SWIcorrYClass
        case 'pval'
            Y_allNets = op.KSpval(SWINetsToInclude);
            SWI_Ylabel = 'pval';
        case 'KSmedian'
            Y_allNets = op.KSMedianDiff(SWINetsToInclude);
            SWI_Ylabel = 'Median shift';
        case 'fracSig'
            Y_allNets = op.fracEventSig(SWINetsToInclude);
            SWI_Ylabel = 'Frac Sig.';
    end
    SWI_Xlabel = 'SWI';
    if analysisParam.corrLogVals
        X_allNets = log10(X_allNets);
        Y_allNets = log10(Y_allNets);
        SWI_Xlabel = [SWI_Xlabel, ' (log10)'];
        SWI_Ylabel = [SWI_Ylabel, ' (log10)'];
    end
    validInds_nets = ~isnan(X_allNets) & ~isinf(X_allNets) & ~isnan(Y_allNets);
    X_allNets = X_allNets(validInds_nets); Y_allNets = Y_allNets(validInds_nets);
    mdl_SWI_netWise = fitlm(squeeze(X_allNets), squeeze(Y_allNets));
    % Calcualte and plot Spearman’s rank correlation coefficient
    [rho,pval] = corr(X_allNets, Y_allNets, 'Type', 'Spearman');
    figure; scatter(X_allNets, Y_allNets, scatterMS) %, 'MarkerStyle', 'o');
    xlabel(SWI_Xlabel); ylabel(SWI_Ylabel)
    title(['p=', num2str(pval, 2), ', \rho=', num2str(rho, 2)], 'FontWeight','Normal');
    hold on; axis manual
    lsl = lsline; lsl.Color='k'; lsl.LineStyle='--';

    % Parameter-wise correlation:
    X_params = squeeze(median(op.SWI, 3));
    X_params = X_params(SWIParamsToInclude);
    switch analysisParam.SWIcorrYClass
        case 'pval'
            Y_params = opComb.KSpval(SWIParamsToInclude);
            SWI_Ylabel = 'pval';
        case 'KSmedian'
            Y_params = opComb.KSMedianDiff(SWIParamsToInclude);
            SWI_Ylabel = 'Median shift';
        case 'fracSig'
            Y_params = opComb.fracEventSig(SWIParamsToInclude);
            SWI_Ylabel = 'Frac Sig.';
    end
    SWI_Xlabel = 'SWI';
    if analysisParam.corrLogVals
        X_params = log10(X_params);
        Y_params = log10(Y_params);
        SWI_Xlabel = [SWI_Xlabel, ' (log10)'];
        SWI_Ylabel = [SWI_Ylabel, ' (log10)'];
    end
    validInds_params = ~isnan(X_params) & ~isinf(X_params) & ~isnan(Y_params);
    X_params = X_params(validInds_params); Y_params = Y_params(validInds_params);
    mdl_SWI_paramWise = fitlm(X_params, Y_params);
    % Calcualte and plot Spearman’s rank correlation coefficient
    [rho,pval] = corr(X_params, Y_params, 'Type', 'Spearman');
    figure; scatter(X_params, Y_params, scatterMS) %, 'MarkerStyle', 'o');
    xlabel(SWI_Xlabel); ylabel(SWI_Ylabel)
    title(['p=', num2str(pval, 2), ', \rho=', num2str(rho, 2)], 'FontWeight','Normal');
    hold on; axis manual
    lsl = lsline; lsl.Color='k'; lsl.LineStyle='--';

    if 0
        % Rank-order linear correlation
        [~, X_sortInds] = sort(X_params);
        [a, X_rank] = sort(X_sortInds);
        [~, Y_sortInds] = sort(Y_params);
        [c, Y_rank] = sort(Y_sortInds);
        mdl_rank = fitlm(X_rank, Y_rank)
        figure; plot(mdl_rank); legend off; title ''
        xlabel('SWI rank'); ylabel('Y rank')
    end
end

myPlotSettings


%% Parameter point-wise SWI correlations

if (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
    Y = op.KSMedianDiff;
    X = op.SWI;
    nEv = op.nEvents;
    op2=struct;
    for ithParam = 1:size(X,1)
        for jthParam = 1:size(X,2)
            %netInds = 1:size(X,3);
            netInds = squeeze(nEv(ithParam,jthParam,:))>=analysisParam.minNEventsSWI;
            if sum(netInds)==0
                op2.r(ithParam,jthParam)=nan;
                op2.p(ithParam,jthParam)=nan;
            else
                [op2.r(ithParam,jthParam),op2.p(ithParam,jthParam)] = corr(squeeze(X(ithParam,jthParam,netInds)), squeeze(Y(ithParam,jthParam,netInds)), 'Rows', 'complete', 'Type', 'Spearman');
            end
        end
    end

    figure
    imagesc(xParamvec, yParamvec, op2.r', 'AlphaData', ~isnan(op2.r')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); %cb.Label.String = 'KSMedianDiff';
    title('SWI correlation r');
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

    figure
    imagesc(xParamvec, yParamvec, op2.p', 'AlphaData', ~isnan(op2.p')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); %cb.Label.String = 'KSMedianDiff';
    title('SWI correlation pval');
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')


    figure
    imagesc(xParamvec, yParamvec, nansum(nEv, 3)', 'AlphaData', ~isnan(nansum(nEv, 3)')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); %cb.Label.String = 'KSMedianDiff';
    title('Total number of events');
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

    figure
    imagesc(xParamvec, yParamvec, sum(nEv>=analysisParam.minNEventsSWI, 3)', 'AlphaData', ~isnan(sum(nEv>=analysisParam.minNEventsSWI, 3)')&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); %cb.Label.String = 'KSMedianDiff';
    title('Number of nets included');
    xlabel(xName,'Interpreter','none'); ylabel(yName,'Interpreter','none')

end

%% op.multiEnvPval

if (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)

    if modelParam.nEnvironments>1
        figure;
        imagesc(xParamvec, yParamvec, opComb.multiEnvPval', 'AlphaData', ~excludedBadConnProb')
        set(gca,'YDir','normal')
        cb = colorbar(); cb.Label.String = 'p-value';
        caxis([0, 0.5])
        xlabel(xName,'Interpreter','none')
        ylabel(yName,'Interpreter','none')
        title('Min(KS test pval) for outlier traj PFs')
    end

    if modelParam.nEnvironments==4
        % KS-test pval comparing Env1 vs Env2 correlation CDFs
        figure;
        imagesc(xParamvec, yParamvec, opComb.envComp_p', 'AlphaData', ~excludedBadConnProb')
        set(gca,'YDir','normal')
        cb = colorbar(); cb.Label.String = 'p-value';
        caxis([0, 0.5])
        xlabel(xName,'Interpreter','none')
        ylabel(yName,'Interpreter','none')
        title('KS-test, Env1 vs Env2')

        % Median shift of correlation CDFs from Env1 minus Env2
        figure;
        imagesc(xParamvec, yParamvec, opComb.envComp_shift', 'AlphaData', ~excludedBadConnProb')
        set(gca,'YDir','normal')
        cb = colorbar(); cb.Label.String = 'median(Env1) - median(Env2)';
        caxis([-0.8*max(abs(opComb.envComp_shift(:))), 0.8*max(abs(opComb.envComp_shift(:)))])
        xlabel(xName,'Interpreter','none')
        ylabel(yName,'Interpreter','none')
        title('Shift in CDF, Env1 vs Env2')
    end
end

%% Functions
function [op, opComb] = loopoOPCalc(op, opComb, analysisParam, simParam, modelParam, PFresultsStruct, resultsStruct)
for ithParam1 = 1:numel(simParam.variedParam(1).range)
    for ithParam2 = 1:numel(simParam.variedParam(2).range)

        % Loop over nets:
        temp.allEnvPvals = [];
        temp.decodeRvals = [];
        temp.allEnvRvals = [];
        temp.allEnvShuffRval = [];
        temp.shuffRvals = [];
        temp.totalnEvents = 0;
        temp.totalnSigEvents = 0;
        for ithNet = 1:simParam.nNets
            [op, temp] = calcOP(op, temp, ithParam1, ithParam2, ithNet, analysisParam, simParam, modelParam, PFresultsStruct, resultsStruct);
        end

        % Cross-env analysis and plots, combined correlations for all networks
        if modelParam.nEnvironments>1
            if analysisParam.compPairWiseEnv
                % Do all all pair-wise CDF comparisons between environment decodes
                pValMat_temp = zeros(size(temp.allEnvRvals, 2), size(temp.allEnvRvals, 2));
                for envI = 1:size(temp.allEnvRvals, 2)
                    for envJ = 1:size(temp.allEnvRvals, 2)
                        [~,p_kstest,~] = kstest2(temp.allEnvRvals(:,envI), temp.allEnvRvals(:,envJ));
                        pValMat_temp(envI, envJ) = p_kstest;
                    end
                end
            else
                % Compare each environment decode CDF to the combined environment CDF
                pValMat_temp = zeros(1, size(temp.allEnvRvals, 2));
                for envI = 1:size(temp.allEnvRvals, 2)
                    [~,p_kstest,~] = kstest2(temp.allEnvRvals(:), temp.allEnvRvals(:,envI));
                    pValMat_temp(envI) = p_kstest;
                end
            end
            if ~isempty(pValMat_temp)
                opComb.multiEnvPval(ithParam1, ithParam2) = min(pValMat_temp(:));
            end
            % Plot for all params, if specified
            if analysisParam.plotCombinedCDFs || isequal(analysisParam.paramPoint, [ithParam1, ithParam2])

                if analysisParam.MSPlots; myPlotSettings(width=1.75, height=1.25, lw=1)
                else; myPlotSettings;
                end
                figure; hold on;
                ecdf(temp.allEnvRvals(:,1)); ecdf(temp.allEnvRvals(:,2)); ecdf(temp.allEnvRvals(:,3)); ecdf(temp.allEnvRvals(:,4));
                [f1,x1]=ecdf(temp.allEnvShuffRval(:,1)); plot(x1,f1,'k')
                [f1,x1]=ecdf(temp.allEnvShuffRval(:,2)); plot(x1,f1,'k')
                [f1,x1]=ecdf(temp.allEnvShuffRval(:,3)); plot(x1,f1,'k')
                [f1,x1]=ecdf(temp.allEnvShuffRval(:,4)); plot(x1,f1,'k')
                %title([simParam.variedParam(1).name, '=', num2str(simParam.variedParam(1).range(ithParam1)), ' ', simParam.variedParam(2).name, '=', num2str(simParam.variedParam(2).range(ithParam2)), ...
                %    ' min(p-val)=', num2str(min(pValMat_temp(:)), 3)])
                legend({'R1', 'L1', 'R2', 'L2', 'Shuffles'}, 'Box', 'off', 'Location', 'Best')
                %legend({'R1', 'L1', 'R2', 'L2', 'R1 shuf.', 'L1 shuf.', 'R2 shuf.', 'L2 shuf.'}, 'Box', 'off', 'Location', 'Best')
                xlabel("|correlation|"); ylabel({'Cumulative', 'proportion'})

                % Plot env correlations of decode correlations
                if (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
                    myPlotSettings
                    figure; plotmatrix(temp.allEnvRvals)
                    title('Correlation of each event decoded by each env')
                end
            end
        end

        % Calculate fraction of events that are significant to one or both
        % environments
        if modelParam.nEnvironments==4 && ~isempty(temp.allEnvPvals)
            significantThreshold = analysisParam.capacitySignificanceAlpha;
            fracSigEnv1 = any(temp.allEnvPvals(:,1:2)<significantThreshold, 2);
            fracSigEnv2 = any(temp.allEnvPvals(:,3:4)<significantThreshold, 2);
            opComb.envFracSigEnv1(ithParam1, ithParam2) = mean(fracSigEnv1);
            opComb.envFracSigEnv2(ithParam1, ithParam2) = mean(fracSigEnv2);
            opComb.envFracSigBoth(ithParam1, ithParam2) = mean(fracSigEnv1 & fracSigEnv2);
            xActual = [1, 1, 2];
            yActual = [mean(fracSigEnv1), mean(fracSigEnv2), mean(fracSigEnv1 | fracSigEnv2)];
            mdl = fitlm(xActual, yActual);
            opComb.envCapacityMdl{ithParam1, ithParam2} = mdl;
            invMdl = fitlm(yActual, xActual);
            minNEvents = 10;
            if numel(fracSigEnv1)>minNEvents
                opComb.envCapacityPred(ithParam1, ithParam2) = predict(invMdl, 1.0);
            end
            if isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) && (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
                if analysisParam.MSPlots
                    myPlotSettings(width=2.25, height=1.5)
                end
                figure;
                hold on
                %plot(mdl)
                x = 0:1:ceil(predict(invMdl, 1.0));
                plot(x, predict(mdl, x'), 'k:')
                scatter(xActual, yActual, 50, 'rx')
                ylim([0, 1])
                yticks([0, 0.5, 1]); yticklabels([0, 50, 100])
                legend off
                %title(['Capacity: param=[', num2str(analysisParam.paramPoint), ']'])
                xlabel('Number of tracks')
                ylabel('Network capacity (%)')
                %xline(predict(invMdl, 1.0))
            end
        else
            disp("Not calculating opComb.envFracSig*")
        end

        % KS-test comparing Env 1 decode vs Env 2 decodes
        if modelParam.nEnvironments==4 &&  ~isempty(temp.allEnvRvals)
            env1Rvals = reshape(temp.allEnvRvals(:,1:2), 1, []);
            env2Rvals = reshape(temp.allEnvRvals(:,3:4), 1, []);
            [~,P] = kstest2(env1Rvals, env2Rvals);
            medianDiff = median(env1Rvals)-median(env2Rvals);
            opComb.envComp_p(ithParam1, ithParam2) = P;
            opComb.envComp_shift(ithParam1, ithParam2) = medianDiff;
        end

        % Combine all net data at current parameter point
        if analysisParam.combineNetData
            if ~isempty(temp.decodeRvals)
                temp.decodeRvals = temp.decodeRvals(randperm(numel(temp.decodeRvals))); % permute, to mix networks' events
                temp.decodeRvals = temp.decodeRvals(1:min(numel(temp.decodeRvals), analysisParam.maxNEvents));
                % numel(temp.decodeRvals)
                [~,p_kstest,KSSTAT] = kstest2(temp.decodeRvals, temp.shuffRvals);
                opComb.KSMedianDiff(ithParam1, ithParam2) = median(temp.decodeRvals)-median(temp.shuffRvals);
                opComb.KSpval(ithParam1, ithParam2) = p_kstest;
                opComb.KSstat(ithParam1, ithParam2) = KSSTAT;
                opComb.fracEventSig(ithParam1, ithParam2) = temp.totalnSigEvents/temp.totalnEvents;
            end
            if (analysisParam.plotCombinedCDFs && ~isempty(temp.decodeRvals)) || isequal(analysisParam.paramPoint, [ithParam1, ithParam2])
                if analysisParam.MSPlots
                    myPlotSettings(width=1.75, height=1.25, lw=1)
                else
                    myPlotSettings;
                end
                figure;
                hold on;
                [f1,x1] = ecdf(temp.shuffRvals);
                [f2,x2] = ecdf(temp.decodeRvals);
                plot(x1, f1, 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
                plot(x2, f2, 'LineWidth', 1, 'color', [0, 0, 0.6])
                %legend({'Shuffles', 'Preplays'}, 'Location', 'Best', 'Box', 'off')
                title(['p=', num2str(p_kstest)], fontsize=10, FontWeight='Normal')
                ylabel({'Cumulative', 'Proportion'})
                xlabel('|correlation|')
            end
        else
            opComb.KSpval(1, ithParam1, ithParam2) = median(op.KSpval(ithParam1, ithParam2, :));
            opComb.KSstat(2, ithParam1, ithParam2) = nanmean(op.KSstat(ithParam1, ithParam2, :));
            opComb.fracEventSig(3, ithParam1, ithParam2) = nanmean(op.totalnSigEvents(ithParam1, ithParam2, :)./op.totalnEvents(ithParam1, ithParam2, :));
        end

    end % ithParam2 loop
end % ithParam1 loop
end

%%

function [op, temp] = calcOP(op, temp, ithParam1, ithParam2, ithNet, analysisParam, simParam, modelParam, PFresultsStruct, resultsStruct)
%
%

resultsExist = ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results);
pMatExists = resultsExist && isfield(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory, 'pMat') && ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat);

if pMatExists

    % Get and format data from resultsStruct
    if analysisParam.useWeightedDecode
        allshuff_r2vals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedR2{:});
        rvals_shuffle = sqrt(abs(vertcat([allshuff_r2vals{:,analysisParam.ithEnv}]'))); % take just analysisParam.ithEnv traj, note: weightedR2 is 1x500, r2 is 500x1
        rvals_allenv_shuffle = sqrt(abs(cell2mat(allshuff_r2vals')'));
        rvals_preplay =        sqrt(abs(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedR2(:,analysisParam.ithEnv)));
        rvals_preplay_allenv = sqrt(abs(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedR2));
        pvals_preplay_allenv = nan(size(rvals_preplay_allenv));
        for ithEvent = 1:size(rvals_preplay_allenv, 1)
            for ithTraj = 1:size(rvals_preplay_allenv, 2)
                pvals_preplay_allenv(ithEvent, ithTraj) = mean( ...
                    rvals_preplay_allenv(ithEvent, ithTraj)<sqrt(abs(allshuff_r2vals{ithEvent, ithTraj})) ...
                    );
            end
        end
        pvals_preplay = pvals_preplay_allenv(:, analysisParam.ithEnv);
        slopes_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedSlope(:,analysisParam.ithEnv);
        slopes_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedSlope{:});
        slopes_shuffle = vertcat(slopes_shuffle{:,analysisParam.ithEnv}); % take just analysisParam.ithEnv traj
    else
        allshuff_r2vals = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_rsquare{:});
        rvals_shuffle = sqrt(abs(vertcat(allshuff_r2vals{:,analysisParam.ithEnv}))); % take just analysisParam.ithEnv traj
        rvals_allenv_shuffle = sqrt(abs(cell2mat(allshuff_r2vals'))');
        rvals_preplay = sqrt(abs(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,analysisParam.ithEnv)));
        rvals_preplay_allenv = sqrt(abs(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare));
        % Equal to frac of shuffle events with greater R2 than actual event's R2
        pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,analysisParam.ithEnv);
        slopes_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.slopes(:,analysisParam.ithEnv);
        try % Struct field added to later simulation
            slopes_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_slopes{:});
            slopes_shuffle = slopes_shuffle(:,analysisParam.ithEnv); % take just analysisParam.ithEnv traj
            slopes_shuffle = cellfun(@transpose,slopes_shuffle,'UniformOutput',false);
            slopes_shuffle = vertcat([slopes_shuffle{:,1}]');
        catch
            slopes_shuffle = nan;
        end
    end

    % Calculate decode variance
    accumMeanDecodeVar = 0;
    for ithEvent = 1:numel(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)
        accumMeanDecodeVar = accumMeanDecodeVar + mean(var(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{analysisParam.ithEnv}.pMat, [], 1));
    end

    % Plot CDFs
    if analysisParam.plotAllCDFs
        figure; hold on; ecdf(rvals_preplay); ecdf(rvals_shuffle); legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
        figure; hold on; ecdf(rvals_preplay_allenv(:,1)); ecdf(rvals_preplay_allenv(:,2)); ecdf(rvals_preplay_allenv(:,3)); ecdf(rvals_preplay_allenv(:,4));
    end

    % KS-test comparing preplay r-vals against shuffles'
    if ~isempty(rvals_preplay) && ~isempty(rvals_shuffle)
        [~,P,KSSTAT] = kstest2(rvals_preplay, rvals_shuffle);
        op.KSpval(ithParam1, ithParam2, ithNet) = P;
        op.KSstat(ithParam1, ithParam2, ithNet) = KSSTAT;
    end

    % Add net's data to struct
    op.KSMedianDiff(ithParam1, ithParam2, ithNet) = median(rvals_preplay) - median(rvals_shuffle);
    op.netGroupID(ithParam1, ithParam2, ithNet) = sub2ind( size(resultsStruct, [1 2]), ithParam1, ithParam2);
    op.decodeSlope(ithParam1, ithParam2, ithNet) = mean(abs(slopes_preplay));
    op.decodeSlopeZScore(ithParam1, ithParam2, ithNet) = [mean(abs(slopes_preplay)) - mean(abs(slopes_shuffle))] / [std(abs(slopes_shuffle))];
    op.decodeEntropy(ithParam1, ithParam2, ithNet) = mean(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.Entropy(:,analysisParam.ithEnv));
    op.decodeVar(ithParam1, ithParam2, ithNet) = accumMeanDecodeVar./numel(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat);
    op.fracEventSig(ithParam1, ithParam2, ithNet) = sum(pvals_preplay<0.05)./numel(pvals_preplay);
    op.nEvents(ithParam1, ithParam2, ithNet) = numel(pvals_preplay);

    % Accumulate all values across nets
    temp.allEnvPvals = [temp.allEnvPvals; pvals_preplay_allenv];
    temp.decodeRvals = [temp.decodeRvals, rvals_preplay'];
    temp.allEnvRvals = [temp.allEnvRvals; rvals_preplay_allenv];
    temp.allEnvShuffRval = [temp.allEnvShuffRval; rvals_allenv_shuffle];
    temp.shuffRvals = [temp.shuffRvals, rvals_shuffle'];
    temp.totalnEvents = temp.totalnEvents + numel(pvals_preplay);
    temp.totalnSigEvents = temp.totalnSigEvents + sum(pvals_preplay<0.05);
end

% Calculate network structure and SWI
% Recreate network
if resultsExist
    E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
    netSeed = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.netSeed;
    netParams=modelParam;
    netParams.spatialBin=2;
    netParams.trackWidth=1;
    netParams.envIDs = 1;% modelParam.envIDs;
    netParams.(simParam.variedParam(1).name) = simParam.variedParam(1).range(ithParam1);
    netParams.(simParam.variedParam(2).name) = simParam.variedParam(2).range(ithParam2);
    netParams = set_depedent_parameters(netParams);
    network = create_network(netParams, 'seed', netSeed);
    if ~all(network.E_indices==E_indices); error('Incorrect network'); end
    W = logical(network.conns(network.E_indices, network.E_indices));
    op.SWI(ithParam1, ithParam2, ithNet) = calcSWI(W);
end

end