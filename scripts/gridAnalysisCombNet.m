%% Plot combined place fields and perform other analyses at a particular parameter point
% gridAnalysisCombNet.m
%
%

%% Choose which simulation results to analyze

decodeFileID = "mainPaperGrid"; % Select which decoding file result ID to analyze

figSettings = "manuscript"; % "manuscript1" for figure R4, "manuscriptAlt" for Figure 7, "manuscript" for other manuscript figures


%% Set up and load data

myPlotSettings

addpath(fullfile("..", "src"))
Config = utils.setConfig;
dataPath = Config.simDataPath;

% Load parameters and results
decodeName = grid.getDecodeFilename(decodeFileID);
analysisParam.loadSpikes = true; % Needs to be true for some extra figures
if analysisParam.loadSpikes
    [modelParam, ...
        simParam, ...
        resultsStruct, ...
        PFresultsStruct, ...
        preplaySpikeTimes_grid, ...
        PFSpikeTimes_grid...
        ] = grid.loadResults(decodeName);
else
    [modelParam, ...
        simParam, ...
        resultsStruct, ...
        PFresultsStruct ...
        ] = grid.loadResults(decodeName);
end


%% Select parameter point(s) to analyze

% analysisParam.paramSetInds = combvec([1:size(resultsStruct, 1)], [1:size(resultsStruct, 2)])'
analysisParam.paramSetInds = combvec([2], [3])'; % best example, for 2022-07-22T18-29 data

analysisParam.paramSetInds_linear = sub2ind([numel(simParam.variedParam(1).range), numel(simParam.variedParam(2).range)], analysisParam.paramSetInds(:,1), analysisParam.paramSetInds(:,2));


%% Set up and analysis parameters
analysisParam.paramPoint = [2, 3]; % Select which example param point to plot
analysisParam.ithEnv = 1; % which environments decoding and place fields to plot/analyze

% Analysis parameters
analysisParam.minPeakRate       = 3;        % Minimum peak PF rate to be considered a place cell
analysisParam.useWeightedDecode = true;     % slope and R2 for correlations by either peak prob or weighted prob
analysisParam.removeBadEvents   = 1;        % Remove decoded events that only had one non-zero probility time bin
analysisParam.nEventsToPlot     = 12;       %
analysisParam.useMeanPFDensity  = true;     %
analysisParam.calcScore         = false;    %

% Plotting options
analysisParam.plotSpikeSequences = false;
analysisParam.plotPresentationFigs = false;
analysisParam.plotExtraPvalMat = false;
analysisParam.plotExtraPlots    = false;	% If true, plot place fields of every network
analysisParam.plotNetsBest      = false;    % For panels 1i-j (plot the best replay event of each network)
analysisParam.plotNetPFs        = false;
analysisParam.plotNetPFsEnv1    = false;    % Plot PFs for the first Env for each network

analysisParam.plotPvalMat       = false;
analysisParam.plotClusterPFs    = false;
analysisParam.plotAllNetStruct  = false;
analysisParam.plotLastNetStruct = false;
analysisParam.plotExtraPFPlots  = false;
analysisParam.plotEventLengthAnalysis = false;
analysisParam.plotPFdistAnalysis = false;
analysisParam.plotMaxJumpCDF	= false;

analysisParam.jumpThreshCm = 5;
analysisParam.nJumpShuffles = 10;

analysisParam.figSettings = figSettings;
% analysisParam.figSettings = 'manuscriptAlt'; Only used for small p-val grid panel
assert(any(strcmp(analysisParam.figSettings, {'standard', 'manuscript', 'manuscript1', 'SfNPoster', 'manuscriptAlt'})))

% PF-dist analysis choices
analysisParam.useEV = 0; % Use mean PF location, rather than PF peak
analysisParam.singleClust = 0; % Only consider cells in single cluster
analysisParam.multipleClust = 0; % Only consider cells in multiple clusters
%analysisParam.endClusts = 1 % only consider cells in the first and last cluster

disp(analysisParam)

% Plot manuscript figures
% Figure 2:
% analysisParam.plotNetsBest = true;
% analysisParam.paramSetInds = combvec([2], [3])';
%
% Figure 4:
%

if analysisParam.plotNetsBest
    assert(analysisParam.loadSpikes)
end


%% Loop over parameter sets
rng('default');
tic
for ithParamSetToPlot = 1:size(analysisParam.paramSetInds, 1)

    ithParamSetToPlot
    ithParam1 = analysisParam.paramSetInds(ithParamSetToPlot,1);
    param1Val = simParam.variedParam(1).range(ithParam1);
    disp(['Param1=', num2str(ithParam1), ': ', simParam.variedParam(1).name, '=', num2str(param1Val)])
    ithParam2 = analysisParam.paramSetInds(ithParamSetToPlot,2);
    param2Val = simParam.variedParam(2).range(ithParam2);
    disp(['Param2=', num2str(ithParam2), ': ', simParam.variedParam(2).name, '=', num2str(param2Val)])

    for i = 1:size(simParam.variedParam, 2)
        modelParam.(simParam.variedParam(i).name) = simParam.variedParam(i).range(analysisParam.paramSetInds(i));
    end


    linearParamInd = find(all([param1Val; param2Val] == simParam.parameterSets_vec, 1));

    if analysisParam.loadSpikes
        spikeTimeCell = preplaySpikeTimes_grid{linearParamInd};
    else
        spikeTimeCell = [];
    end
    op = calcOP(analysisParam, simParam, modelParam, resultsStruct, PFresultsStruct, spikeTimeCell, ithParam1, ithParam2);


    %% Plot histogram of place field map correlations (all comparisons for all nets)
    if analysisParam.plotNetPFs
        switch analysisParam.figSettings
            case 'standard'
                myPlotSettings
            case 'manuscript'
                myPlotSettings(width=2.5, height=1.5)
            case 'SfNPoster'
                myPlotSettings(width=3, height=2.5, lw=3, afs=24, alw=3)
        end

        withinEnvCorr = [op.allNetPFEnvCorrs(:,2,1), op.allNetPFEnvCorrs(:,4,3)];
        acrossEnvCorr =  op.allNetPFEnvCorrs(:,3:4,1:2);
        allCorr = [withinEnvCorr(:); acrossEnvCorr(:)];
        %figure; histogram(op.allNetPFEnvCorrs); ylabel('All pair-wise comparisons (count)'); xlabel('PF map correlation')
        %figure; histogram(withinEnvCorr); ylabel('Within-env (count)'); xlabel('PF map correlation')
        %figure; histogram(acrossEnvCorr); ylabel('Across-env (count)'); xlabel('PF map correlation')

        nBins = 2*ceil(sqrt(numel([withinEnvCorr(:); acrossEnvCorr(:)])));
        nBins = 15;
        % figure; histogram(allCorr, nBins); ylabel('All pair-wise comp. (count)'); xlabel('PF map correlation')

        [N,EDGES] = histcounts(allCorr, nBins);
        figure; hold on
        histogram(withinEnvCorr, EDGES)
        histogram(acrossEnvCorr, EDGES)
        legend({'Within-env.','Across-env.'}, 'location', 'north', 'Box', 'off')
        ylabel('All pair-wise (count)'); xlabel('Remapping correlation (r)')
        %xlim([-inf, 1])

    end

    %% Misc. extra plots
    if analysisParam.plotExtraPlots
        % timebins sum distribution
        figure; histogram( round(op.allDecodeSums, 10), 'Normalization', 'Probability')
        xlabel('Prob. sum in time-bin'); ylabel('Prob. (all decode time-bins)')

        % Decode correlation is negatively correlated with p-value
        figure; scatter(op.allEventR2s, op.allEventPvals); xlabel('regression R'); ylabel('p-val against shuffles')

        % Plot mean decode position
        figure; plot(op.avgDecodePos./sum(op.avgDecodePos))
        xlabel('Position (2 cm bins)'); ylabel('Probability Density')
    end

    %% Plot ECDF of maxJump values
    if analysisParam.plotMaxJumpCDF
        switch analysisParam.figSettings
            case 'standard'
                myPlotSettings
            case 'manuscript'
                myPlotSettings(width=2.5,height=1.5)
            case 'SfNPoster'
                % myPlotSettings(width=4, height=2.75, lw=1.5, fzs=14, alw=1.25) % ppt format
                % myPlotSettings(width=3, height=1.5) % for poster
                myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3) % SfN-poster format
        end
        if ~isempty(op.allEventMaxJumps)
            allEventMaxJumps_ecdf = (op.allEventMaxJumps);
            maxJump_shuff_ecdf= (vertcat([op.allShuffleMaxJumps{:}]'));
            warning('Need to fix JD ecdf plots')
            %figure; hold on; ecdf(sqrt(allEventMaxJumps_ecdf)); ecdf(real(sqrt(maxJump_shuff_ecdf))); %ecdf(sqrt(r2vals_shuff_ecdf));
            figure; hold on; cdfplot(sqrt(allEventMaxJumps_ecdf)); cdfplot(real(sqrt(maxJump_shuff_ecdf))); %ecdf(sqrt(r2vals_shuff_ecdf));
            legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
            xlabel('Max Jump (proportion of track)'); ylabel({'Cumulative','proportion'});
            [H,P,KSSTAT] = kstest2(sqrt(abs(allEventMaxJumps_ecdf)), sqrt(abs(maxJump_shuff_ecdf)) )
            title([simParam.variedParam(1).name, '=', num2str(simParam.variedParam(1).range(ithParam1)), ' ', ...
                simParam.variedParam(2).name, '=', num2str(simParam.variedParam(2).range(ithParam2)), ...
                ' pval=', num2str(P), ...
                ' nEvents=', num2str(numel(allEventMaxJumps_ecdf))])

            title ''
            legend({'Preplay', 'Shuffle'}, 'Location', 'Best'); legend boxoff
            grid off

            h = get(gca,'children');
            %    set(h(1), 'LineWidth', 1, 'color', [0.7, 0, 0], 'LineStyle', '-.')
            set(h(1), 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
            set(h(2), 'LineWidth', 1, 'color', [0, 0, 0.6])
            %set(h(2), 'LineWidth', 1, 'color', 'k', 'LineStyle', '--'); set(h(3), 'LineWidth', 1, 'color', [0, 0, 0.6])
        else
            disp('No events to plot jump CDF')
        end
    end


    %% Plot p-value matrix
    if analysisParam.plotPvalMat
        switch analysisParam.figSettings
            case 'standard'
                myPlotSettings
            case 'manuscript'
                myPlotSettings(width=2.5,height=1.5)
            case 'manuscript1'
                myPlotSettings(width=1.5,height=1)
            case 'manuscriptAlt'
                myPlotSettings(width=1.8, height=1.25)
            case 'SfNPoster'
                myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3) % SfN-poster format
        end

        disp('Starting pvalMat analysis...')

        eventRs     = sqrt(op.allEventR2s);
        shuffleRs   = cellfun(@sqrt, op.allShuffleR2s, 'UniformOutput', false);
        eventJD     = op.allEventMaxJumps;
        shuffleJD   = op.allShuffleMaxJumps;
        [figHandle, pValMat] = utils.makePvalMatrix(eventRs, shuffleRs, ...
            eventJD, shuffleJD, plotExtra=analysisParam.plotExtraPvalMat, nanThreshIs0=true);

    end


    %% Plot ECDF of r values
    switch analysisParam.figSettings
        case 'standard'
            myPlotSettings
        case 'manuscript'
            myPlotSettings(width=2.5,height=1.5)
        case 'manuscript1'
            myPlotSettings(width=1,height=1)
        case 'SfNPoster'
            myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3) % SfN-poster format
    end

    if ~isempty(op.allEventR2s)
        allEventR2s_ecdf = (op.allEventR2s);
        r2vals_shuff_ecdf= (vertcat(op.allShuffleR2s{:}));

        %figure; hold on; ecdf(sqrt(allEventR2s_ecdf)); ecdf(real(sqrt(r2vals_shuff_ecdf))); %ecdf(sqrt(r2vals_shuff_ecdf));
        figure; hold on; cdfplot(sqrt(allEventR2s_ecdf)); cdfplot(real(sqrt(r2vals_shuff_ecdf))); %ecdf(sqrt(r2vals_shuff_ecdf));
        % figure; hold on; ecdf(allEventR2s_ecdf); ecdf(rvals_ecdf)
        legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
        xlabel('|Correlation|'); ylabel({'Cumulative','proportion'});
        [~,P,KSSTAT] = kstest2(sqrt(abs(allEventR2s_ecdf)), sqrt(abs(r2vals_shuff_ecdf)) );
        disp(['Decode pval: ', num2str(P)])
        disp(['Decode ks-stat: ', num2str(KSSTAT)])
        title([simParam.variedParam(1).name, '=', num2str(simParam.variedParam(1).range(ithParam1)), ' ', ...
            simParam.variedParam(2).name, '=', num2str(simParam.variedParam(2).range(ithParam2)), ...
            ' pval=', num2str(P), ...
            ' nEvents=', num2str(numel(allEventR2s_ecdf))])
        title ''
        legend({'Preplay', 'Shuffle'}, 'Location', 'Best'); legend boxoff
        grid off

        h = get(gca,'children');
        %    set(h(1), 'LineWidth', 1, 'color', [0.7, 0, 0], 'LineStyle', '-.')
        set(h(1), 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
        set(h(2), 'LineWidth', 1, 'color', [0, 0, 0.6])
        %set(h(2), 'LineWidth', 1, 'color', 'k', 'LineStyle', '--'); set(h(3), 'LineWidth', 1, 'color', [0, 0, 0.6])
    else
        disp('No events to plot r-val CDF')
    end


    %% Plot best sequences
    switch analysisParam.figSettings
        case 'standard'
            myPlotSettings(width=6, height=4)
        case 'manuscript'
            myPlotSettings(width=3.5, height=1.75, ttlfsz=1.0)
        case 'SfNPoster'
            myPlotSettings(width=7, height=7, lw=3, fzs=24, alw=3) % SfN-poster format, Too big for location for secondary
    end

    % Decodes
    tBinSz = resultsStruct(ithParam1, ithParam2, 1).results.replaytrajectory.tBinSz;
    if ~all(isinf(op.bestEventpVals))
        if analysisParam.nEventsToPlot>0
            [pvals_sorted, eventSortInd] = sort(op.bestEventpVals);
            figure; tfig = tiledlayout('flow', 'Padding', 'none', 'TileSpacing', 'compact');
            title(tfig, ['Best ', num2str(numel(op.bestEvents_decode)), ' of ', num2str(op.nEvents), ' events'])
            for ithEventToPlot = eventSortInd'
                nexttile

                [yt, xt] = size(op.bestEvents_decode{ithEventToPlot});
                imagesc([1:xt]*(tBinSz), [1:yt]*(modelParam.spatialBin*100), op.bestEvents_decode{ithEventToPlot})
                colormap hot
                title(['pval=', num2str(op.bestEventpVals(ithEventToPlot))], 'FontWeight', 'normal')
                title(['r=', num2str(op.bestEventAbsrVals(ithEventToPlot), 2), '; jd=', num2str(op.bestEventjdVals(ithEventToPlot), 2)], 'FontWeight', 'normal')
                % caxis(([0, 0.75*max(op.bestEvents_decode{ithEventToPlot}, [], 'all')]))
                caxis(([0, 0.25]))

                set(gca,'ytick',[])
                %yt = yticks; yticklabels(yt*(modelParam.spatialBin*100)); ylabel('Position (cm)')
                %xt = xticks; xticklabels(xt*(tBinSz)); xlabel('Time (ms)')
                %xticks('auto')
                %if find(ithEventToPlot==eventSortInd)<=(2*analysisParam.nEventsToPlot/3)
                %    set(gca,'xtick',[])
                %end
                set(gca,'xtick',[])
            end
        end
        title(tfig, '')

        if ~exist('bestDecodeColorBarFigHandle', 'var')
            bestDecodeColorBarFigHandle = figure;
            title('colorbar', 'fontSize', 19); colorbar; colormap(hot); caxis([0, 0.25])
        end

        % Relative ranks (corresponding to decodes)
        if analysisParam.nEventsToPlot>0 && analysisParam.plotSpikeSequences
            [pvals_sorted, eventSortInd] = sort(op.bestEventpVals);
            figure; tfig = tiledlayout('flow');
            %title(tfig, ['Best ', num2str(numel(op.bestEvents_relRank)), ' of ', num2str(op.nEvents), ' events'])
            disp(['Best ', num2str(numel(op.bestEvents_relRank)), ' of ', num2str(op.nEvents), ' events'])
            for ithEventToPlot = eventSortInd'
                nexttile
                eventSeq = op.bestEvents_relRank{ithEventToPlot};
                scatter(1:sum(~isnan(eventSeq)), eventSeq(~isnan(eventSeq)), [], op.bestEvents_relRankColor{ithEventToPlot}(~isnan(eventSeq)), 'filled')
                %scatter(1:modelParam.n_E, op.bestEvents_relRank{ithEventToPlot})
                title(['pval=', num2str(op.bestEventpVals(ithEventToPlot))])
                %caxis(([0, 0.5*max(bestEvents{ithEventToPlot}, [], 'all')]))
            end
        end

        % Spike raster, corresponding to above Relative Ranks plot
        if analysisParam.nEventsToPlot>0 && ~isempty(op.bestEvents_raster{1}) && analysisParam.plotSpikeSequences
            MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
            [pvals_sorted, eventSortInd] = sort(op.bestEventpVals);
            figure; tfig = tiledlayout('flow');
            %title(tfig, ['Best ', num2str(numel(op.bestEvents_relRank)), ' of ', num2str(op.nEvents), ' events'])
            for ithEventToPlot = eventSortInd'
                nexttile
                eventRaster = op.bestEvents_raster{ithEventToPlot};
                plotColor =  []; % op.bestEvents_relRankColor{ithEventToPlot}(~isnan(eventSeq));
                % imagesc(eventRaster)
                if ~isempty(eventRaster)
                    % figure; plotSpikeRaster( netPreplayRaster, 'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat);
                    plotSpikeRaster( eventRaster, ...
                        'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
                        'MarkerFormat', MarkerFormat);
                    set(gca,'YDir','normal')
                    title(['pval=', num2str(op.bestEventpVals(ithEventToPlot))])
                end
            end
        end

    else
        disp('No events to plot best examples')
    end


    %% Plot example shuffled events

    if analysisParam.plotPresentationFigs
        myPlotSettings(width=1.5, height=1.25)
        eventToPlot = 11;
        eventMat = op.bestEvents_decode{eventToPlot};
        [yt, xt] = size(eventMat);
        figure
        imagesc([1:xt]*(tBinSz), [1:yt]*(modelParam.spatialBin*100), eventMat)
        colormap hot
        caxis(([0, 0.25]))
        ylim([0, 100])
        %xlim([5, 50])
        set(gca,'ytick',[0, 100])
        set(gca,'xtick',[10, 50])
        xlabel('Time')
        ylabel('Position (cm)')

        myPlotSettings(width=3, height=2.5)
        rng(123)
        figure; tiledlayout('flow', TileSpacing='compact', padding='compact'); % loose, compact, tight, none
        nShuf = 9;
        for ithShuf = 1:nShuf
            nexttile
            imagesc([1:xt]*(tBinSz), [1:yt]*(modelParam.spatialBin*100), eventMat(:,randsample(1:xt,xt)))
            colormap hot
            caxis(([0, 0.25]))
            set(gca,'ytick',[])
            set(gca,'xtick',[])
        end
    end


    %% Plot combined place fields
    switch analysisParam.figSettings
        case 'standard'
            myPlotSettings(width=6, height=4)
        case 'manuscript'
            myPlotSettings(width=1.75, height=1.5)
        case 'SfNPoster'
            % myPlotSettings(width=4, height=3, lw=2, afs=14, alw=2) % ppt format
            myPlotSettings(8, 7, 3, 24, [], [], 3) % SfN-poster format
    end

    row_all_zeros1 = find(all( op.PFmatE_all==0, 2)) ;
    row_n_all_zeros1 = find(~all( op.PFmatE_all==0, 2)) ;
    [peakRate,peakRateLocation] = max(squeeze(op.PFmatE_all(row_n_all_zeros1,:)), [], 2);
    [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    PFpeaksSequence = [row_n_all_zeros1(sortedCellIndsbyPeakRateLocation); row_all_zeros1];
    [peakRate, peakRateLocation_all] = max(op.PFmatE_all, [], 2);

    normRates = 1;
    if normRates
        rateDenomAll = max(op.PFmatE_all(PFpeaksSequence,:), [], 2);
        caxmax = 1;
    else
        rateDenomAll = 1;
        caxmax = max(op.PFmatE_all, [], 'all');
    end

    PFdatatoplot = op.PFmatE_all(PFpeaksSequence(rateDenomAll>analysisParam.minPeakRate),:)./rateDenomAll(rateDenomAll>analysisParam.minPeakRate);
    figure;
    imagesc( fliplr(PFdatatoplot(:,5:end)) );
    title('All nets');
    colorbar; caxis([0, caxmax])
    xlabel('Position (2 cm bin)');
    ylabel('Cell (sorted)');
    colormap(jet); title ''; colorbar off
    xt = xticks; xticklabels(xt*(modelParam.spatialBin*100)); xlabel('Position (cm)')
    % set(gca,'xtick',[]); xlabel('')

    if analysisParam.plotPresentationFigs
        % Plot shuffled cell-identity place fields
        shuffleInds = randsample(size(PFdatatoplot, 1), size(PFdatatoplot, 1));
        figure;
        imagesc( fliplr(PFdatatoplot(shuffleInds,5:end)) );
        title('All nets');
        colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)');
        ylabel('Cell (shuffled)');
        colormap(jet); title ''; colorbar off
        xt = xticks; xticklabels(xt*(modelParam.spatialBin*100)); xlabel('Position (cm)')
    end

    %% Plot cluster-wise place fields
    if analysisParam.plotClusterPFs
        figure; tiledlayout(modelParam.clusters,1);
        singularMembership = 1;
        for ithCluster = modelParam.clusters:-1:1
            if singularMembership
                isClusterMember = [sum(op.cluster_matAll_sorted(PFpeaksSequence(rateDenomAll>analysisParam.minPeakRate), :), 2)==1] .* [op.cluster_matAll_sorted(PFpeaksSequence(rateDenomAll>analysisParam.minPeakRate),ithCluster)==1];
            else
                isClusterMember = [op.cluster_matAll_sorted(PFpeaksSequence(rateDenomAll>analysisParam.minPeakRate),ithCluster)==1];
            end
            clusterPF = isClusterMember .* op.PFmatE_all(PFpeaksSequence(rateDenomAll>analysisParam.minPeakRate),:)./rateDenomAll(rateDenomAll>analysisParam.minPeakRate);
            zeroInds = all(clusterPF==0, 2);
            nexttile; imagesc( clusterPF(~zeroInds,:) );
            %xaxis('off')
            set(gca,'xtick',[])

            % title('All nets'); colorbar; caxis([0, caxmax])
            % xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        end
    end


    %% Extra plots
    if analysisParam.plotExtraPFPlots
        myPlotSettings
        % Example network PFs
        for ithNet = 1%:10
            cellsToPlot = [rateDenomAll>analysisParam.minPeakRate] & [op.netInd_all(PFpeaksSequence)==ithNet];
            figure; imagesc( op.PFmatE_all(PFpeaksSequence(cellsToPlot),:)./rateDenomAll(cellsToPlot)); title('Example net'); colorbar; caxis([0, caxmax])
            xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        end

        % myPlotSettings(height=3, width=1.5) % for poster
        % myPlotSettings(height=4, width=3.5, lw=3, afs=24, alw=3) % SfN-poster format

        % 'Peak Rate'
        peakRate= max(op.PFmatE_all, [], 2);
        figure; histogram(peakRate(peakRate>analysisParam.minPeakRate), 20); xlabel('Place field peak (Hz)'); ylabel('Place cells (count)');
        box off

        % 'kstest'
        ksstat_vec = [];
        for i = 1:size(op.PFmatE_all, 1)
            [~,p,ksstat,~] = kstest( ( op.PFmatE_all(i,:)-mean(op.PFmatE_all(i,:), 2) )./(std(op.PFmatE_all(i,:), [], 2)+eps  ) );
            ksstat_vec(i) = ksstat;
        end
        figure; histogram(ksstat_vec(peakRate>analysisParam.minPeakRate)); xlabel('kstest stat'); ylabel('Place cells (count)');

        % 'sparsity'
        cellSparsity =  mean( op.PFmatE_all>[0.25*max(op.PFmatE_all, [], 2)], 2 );
        figure; histogram(cellSparsity(peakRate>analysisParam.minPeakRate)); xlabel('PF sparsity'); ylabel('Place cells (count)');

        % 'information'
        spatialInfo = nanmean( [op.PFmatE_all./mean(op.PFmatE_all, 2)] .* log(( op.PFmatE_all+eps )./mean(op.PFmatE_all, 2) ), 2 );
        figure; histogram(spatialInfo(peakRate>analysisParam.minPeakRate), 20); xlabel('Information (bits/s)'); ylabel('Place cells (count)');
        box off


        figure; scatter(cellSparsity, spatialInfo); xlabel('PF sparsity'); ylabel('Spatial info.')
        figure; scatter(peakRate, spatialInfo); xlabel('Peak rate (Hz)'); ylabel('Spatial info.')

        mdl1 = fitlm(op.nClustMemb_all, spatialInfo); figure; plot(mdl1); xlabel('n cluster membership'); ylabel('Spatial information');
        title(['pval=', num2str(mdl1.Coefficients.pValue(2), 2), ', rsqr=', num2str(mdl1.Rsquared.Ordinary, 2)]); legend off; mdl1
        mdl2 = fitlm(op.nClustMemb_all, cellSparsity); figure; plot(mdl2); xlabel('n cluster membership'); ylabel('Spatial sparsity');
        title(['pval=', num2str(mdl2.Coefficients.pValue(2), 2), ', rsqr=', num2str(mdl2.Rsquared.Ordinary, 2)]); legend off; mdl2
        mdl3 = fitlm(op.nClustMemb_all, peakRate); figure; plot(mdl3); xlabel('n cluster membership'); ylabel('Peak rate (Hz)');
        title(['pval=', num2str(mdl3.Coefficients.pValue(2), 2), ', rsqr=', num2str(mdl3.Rsquared.Ordinary, 2)]); legend off; mdl3

        cellsToPlot = [spatialInfo>1.87]; % [rateDenomAll>2];
        figure; imagesc( op.PFmatE_all(PFpeaksSequence(cellsToPlot),:)./rateDenomAll(cellsToPlot)); title('All nets'); colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');


        % 'Example best place fields'
        [~, SIind] = sort( ksstat_vec' .* [peakRateLocation_all>10&peakRateLocation_all<40] .* [peakRate>4] .* [cellSparsity<0.5], 'descend' );
        % [~, SIind] = sort( spatialInfo' .* [peakRateLocation_all>10&peakRateLocation_all<40] .* [meanPeakRate>4] .* [cellSparsity>0.5], 'descend' );
        figure; plot(1:size(op.PFmatE_all, 2), op.PFmatE_all(SIind(1:4),:))
        xlabel('Location (2 cm bins)'); ylabel('Firing rate (Hz)'); title('Good example place fields')
        figure; plot(1:size(op.PFmatE_all, 2), op.PFmatE_all(SIind(5:8),:))
        xlabel('Location (2 cm bins)'); ylabel('Firing rate (Hz)'); title('Good example place fields')
        figure; plot(1:size(op.PFmatE_all, 2), op.PFmatE_all(SIind(9:12),:))
        xlabel('Location (2 cm bins)'); ylabel('Firing rate (Hz)'); title('Good example place fields')


        % Mean Place field activity
        figure; plot(mean(fliplr(PFdatatoplot), 1))

        % Distribution of PF peaks
        figure; plot(mean(fliplr(PFdatatoplot==1)))
    end


    %% Plot PF score, if it was calculated
    if analysisParam.calcScore
        figure; histogram(op.PFscores_all); xlabel('Place field score (lower is better)'); ylabel('Network (count)');
    end


    %% Correlation of event length with decode correlation

    if analysisParam.plotEventLengthAnalysis
        myPlotSettings
        mdl_preplay = fitlm(op.allEventLengths, op.allEventR2s); figure; plot(mdl_preplay); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)');
        title(['pval=', num2str(mdl_preplay.Coefficients.pValue(2), 2), ', rsqr=-', num2str(sqrt(mdl_preplay.Rsquared.Ordinary), 2)]); legend off; mdl_preplay
        figure; scatterhist(op.allEventLengths, op.allEventR2s); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)');

        mdl_shuff = fitlm(op.allShuffleLengths, vertcat(op.allShuffleR2s{:}) ); figure; plot(mdl_shuff); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)');
        title(['pval=', num2str(mdl_shuff.Coefficients.pValue(2), 2), ', rsqr=-', num2str(sqrt(mdl_shuff.Rsquared.Ordinary), 2)]); legend off; mdl_shuff
        figure; scatterhist(op.allShuffleLengths, vertcat(op.allShuffleR2s{:})); xlabel('Event time bins (count)'); ylabel('Decode correlation (|r|)');
    end


    %% Plot PF-distance dependent connection probability and analysis

    if analysisParam.plotPFdistAnalysis
        myPlotSettings
        xPosVals = [0:(numel(net_rand_posDepConProb)-1)];

        figure; hold on; title(['All networks'])
        plot(xPosVals, mean(op.allPosDepConProb_norm, 1), 'r');
        xlabel('PF peak distance (2 cm bins)'); ylabel('P(connection | PF dist.)');

        mdl = fitlm(xPosVals, mean(op.allPosDepConProb_norm, 1));
        figure; plot(mdl)

        for ithNet = 1:3%size(resultsStruct, 3)
            figure; hold on; title(['Network #', num2str(ithNet)])
            plot(xPosVals, op.allPosDepConProb_norm(ithNet,:), 'r');
            xlabel('PF peak distance (2 cm bins)'); ylabel('P(connection | PF dist.)');
        end

        %{
        % Plot against shuffles
        figure; hold on; title(['All networks'])
        plot(xPosVals, mean(op.allRand_posDepConProb, 1), 'k:'); plot(xPosVals, mean(op.allRand2_posDepConProb, 1), 'b--');
        plot(xPosVals, mean(op.allPosDepConProb, 1), 'r');
        legend({'Shuffle peak-dist.', 'Shuffle conn.-prob.', 'Actual'}, 'Location', 'Best')
        xlabel('PF peak distance (2 cm bins)'); ylabel('Connection prob.');
        
        for ithNet = 1:3%size(resultsStruct, 3)
            figure; hold on; title(['Network #', num2str(ithNet)])
            plot(xPosVals, op.allRand_posDepConProb(ithNet,:), 'k:'); plot(xPosVals, op.allRand2_posDepConProb(ithNet,:), 'b--');
            plot(xPosVals, op.allPosDepConProb(ithNet,:), 'r');
            legend({'Shuffle peak-dist.', 'Shuffle conn.-prob.', 'Actual'}, 'Location', 'Best')
            xlabel('PF peak distance (2 cm bins)'); ylabel('Connection prob.');
        end

        % Compare shuffle and nonshuffed
          X = [op.allPosDepConProb; op.allRand2_posDepConProb; op.allRand_posDepConProb]';
        % X = [op.allPosDepConProb; op.allRand2_posDepConProb]';
        [p,tbl,stats] = anova2(X, size(op.allPosDepConProb, 1));
        figure; [c,m,h,gnames] = multcompare(stats)
        %}

        %{
        % Calculate correlation for each net individually
        for ithNet = 1:10
            mdl = fitlm(xPosVals, op.allPosDepConProb_norm(ithNet, :));
            mdl.Coefficients.Estimate(2)
            mdl.Coefficients.pValue(2)
        end
        %}

        %{
        % Calculate ANOVA2 and post-hoc
        X = [op.allPosDepConProb; op.allRand2_posDepConProb];
        [p,tbl,stats] = anova2(X, size(op.allPosDepConProb, 1))
        figure; [c,m,h,gnames] = multcompare(stats)
        %}
    end

    %% Plot number of decode jumps

    switch analysisParam.figSettings
        case 'standard'
            myPlotSettings
        case 'manuscript'
            myPlotSettings(width=3, height=2, lw=2)
        case 'SfNPoster'
            myPlotSettings(width=3, height=2.5, lw=3, afs=24, alw=3)
    end

    [H, p_kstest, KSSTAT] = kstest2(op.nJumps, op.nJumpsShuffled);
    % figure; hold on; ecdf(op.nJumps); ecdf(op.nJumpsShuffled);
    figure;
    hold on;
    [f1,x1] = ecdf(op.nJumpsShuffled);
    [f2,x2] = ecdf(op.nJumps);
    plot(x1, f1, 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
    plot(x2, f2, 'LineWidth', 1, 'color', [0, 0, 0.6])
    title(['p=', num2str(p_kstest), ', jump thresh (cm): ', num2str(analysisParam.jumpThreshCm)] ...
        , fontsize=10, FontWeight='Normal')
    ylabel({'Cumulative', 'Proportion'})
    xlabel('Number of steps')
    xline(mean(op.nJumpsShuffled), 'color', [1, 0.3, 0.3], 'LineStyle', '-.');
    xline(mean(op.nJumps), 'color', [0, 0, 0.6]);
    legend({'Shuffles', 'Preplays', num2str(mean(op.nJumpsShuffled)), num2str(mean(op.nJumps))}, ...
        'Location', 'Best', 'Box', 'off')

end

runTime = toc;
disp(['Runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS')])



%% Functions

function op = calcOP(analysisParam, simParam, modelParam, resultsStruct, PFresultsStruct, spikeTimeCell, ithParam1, ithParam2)
% calcOP
%
% Runs the analysis across networks at the specified parameter point and
% returns the results in the output struct op.
%

ithTrial = 1; % Some of the analyses might assume this is true

op.PFmatE_all = [];
op.netInd_all = []; % Index of which net each cell came from
op.nClustMemb_all = []; % number of clusters each cell is a member of
op.PFscores_all = [];
op.bestEventpVals = inf(analysisParam.nEventsToPlot, 1); % analysisParam.nEventsToPlot best event pvals
op.bestEventAbsrVals = inf(analysisParam.nEventsToPlot, 1); % analysisParam.nEventsToPlot best event pvals
op.bestEventjdVals = inf(analysisParam.nEventsToPlot, 1); % analysisParam.nEventsToPlot best event pvals
op.bestEvents_decode = cell(analysisParam.nEventsToPlot, 1);
op.bestEvents_relRank = cell(analysisParam.nEventsToPlot, 1);
op.bestEvents_relRankColor = cell(analysisParam.nEventsToPlot, 1);
op.bestEvents_raster = cell(analysisParam.nEventsToPlot, 1);
op.nEvents = 0;
op.cluster_matAll_sorted = [];
op.allDecodeSums = []; % For all time bins of all events of all nets, what is the sum of the bins probability
op.allEventR2s = [];
op.allEventPvals = [];
op.allEventMaxJumps = [];
op.allShuffleR2s = [];
op.allShuffleMaxJumps = [];
op.allEventLengths = [];
op.allShuffleLengths = [];
op.allPosDepConProb = zeros(size(resultsStruct, 3),     size(modelParam.gridxvals, 2)); % Position depenedent connection probability
op.allRand_posDepConProb = zeros(size(resultsStruct, 3), size(modelParam.gridxvals, 2)); % Rand-pos Position depenedent connection probability
op.allRand2_posDepConProb = zeros(size(resultsStruct, 3), size(modelParam.gridxvals, 2)); % Rand-pos Position depenedent connection probability
op.allPosDepConProb_norm = zeros(size(resultsStruct, 3), size(modelParam.gridxvals, 2)); % Position depenedent connection probability
op.avgDecodePos = zeros(size(modelParam.gridxvals));
op.nJumps = [];
op.nJumpsShuffled = [];
op.allNetPFEnvCorrs = nan(modelParam.nNets, modelParam.nEnvironments, modelParam.nEnvironments);

for ithNet = 1:size(resultsStruct, 3)

    % Get matrix of PFs
    if iscell(PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields)
        % if modelParam.nEnvironments>1
        PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields{analysisParam.ithEnv};
    elseif ismatrix(PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields) && analysisParam.ithEnv==1
        %else
        PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;
    else
        disp(['PF data error ', num2str(ithParam1), ' ', num2str(ithParam2), ' ', num2str(ithNet)])
        PFmat = [];
        E_indices = [];
    end
    E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
    netSeed = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.netSeed;

    op.PFmatE_all = [op.PFmatE_all; PFmat(E_indices,:)];
    op.netInd_all = [op.netInd_all; repmat(ithNet, numel(E_indices), 1)];

    % Recreate network
    netParams=modelParam;
    netParams.spatialBin=2;
    netParams.trackWidth=2;
    netParams.envIDs = modelParam.envIDs;% modelParam.envIDs;
    for i = 1:size(simParam.variedParam, 2)
        netParams.(simParam.variedParam(i).name) = simParam.variedParam(i).range(analysisParam.paramSetInds(i));
    end
    netParams = set_depedent_parameters(netParams);
    network = create_network(netParams, 'seed', netSeed);
    if ~all(network.E_indices==E_indices); error('Incorrect network'); end
    if mod(analysisParam.ithEnv, 2)==1
        [~, clusterOrder] = sort(network.clusterSequence(:,analysisParam.ithEnv));
    else
        [~, clusterOrder] = sort( [netParams.clusters+1] - flipud(network.clusterSequence(:,analysisParam.ithEnv)) );
    end
    sortedClusters = network.cluster_mat(clusterOrder,E_indices)';

    op.nClustMemb_all = [op.nClustMemb_all; sum(network.cluster_mat(:,E_indices), 1)'];
    op.cluster_matAll_sorted = [op.cluster_matAll_sorted; sortedClusters] ;


    % Load trial spikes
    if analysisParam.loadSpikes
        netPreplaySpikeTimes = squeeze(spikeTimeCell(ithNet,ithTrial,:));
        netPreplayRaster = zeros(modelParam.n, modelParam.t_steps_preplay);
        for ithCell = 1:modelParam.n
            cellSpikeInds = round(netPreplaySpikeTimes{ithCell}/modelParam.dt);
            netPreplayRaster(ithCell,cellSpikeInds) = 1;
        end
        netPreplayRaster = logical(netPreplayRaster);
        % Plot whole-trial raster
        if 0
            MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
            % figure; plotSpikeRaster( netPreplayRaster, 'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat);
            figure; plotSpikeRaster( [netPreplayRaster(network.E_indices,:); ...
                netPreplayRaster(network.I_indices,:) ], ...
                'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
                'MarkerFormat', MarkerFormat);
            yline(numel(network.E_indices))
        end
    end


    %% Calculate and plot all pair-wise PF correlations

    if isequal(ithNet, 1)
        myPlotSettings(width=2, heigh=1.25)
        pfMatE = PFmat(E_indices,:);
        clusterMatE = network.cluster_mat(:,E_indices);
        nECells = size(pfMatE, 1);
        opSame = [];
        opDiff = [];
        for i = 1:(nECells-1)
            for j = i+1:nECells
                pairWisecorr = corr(pfMatE(i,:)', pfMatE(j,:)');
                if any(clusterMatE(:,i) & clusterMatE(:,j))
                    opSame = [opSame, pairWisecorr];
                else
                    opDiff = [opDiff, pairWisecorr];
                end
            end
        end
        % figure; hold on; plot(pfMatE(1,:)); plot(pfMatE(2,:))
        figure;
        hold on;
        [f1, x1] = ecdf(opSame);
        plot(x1, f1)
        if ~isempty(opDiff)
            [f2, x2] = ecdf(opDiff);
            plot(x2, f2)
        end
        xlabel('Pairwise PF correlation'); ylabel({'Cumulative', 'proportion'})
        %legend({'Shared cluster', 'No shared cluster'}, 'Location', 'Best')
        %figure; histogram(opSame)
        %figure; histogram(opDiff)
        figure; hold on; h = zeros(2, 1);
        h(1) = plot(NaN,NaN); h(2) = plot(NaN,NaN);
        legend(h, {'Shared cluster', 'No shared cluster'}, 'Location', 'Best');
        axis off
    end


    %% Calculate all pair-wise place field distances and their connectivity
    nValidPairs = 0;
    net_posDepConProb = zeros(size(modelParam.gridxvals)); % Position depenedent connection probability
    net_rand_posDepConProb = zeros(size(modelParam.gridxvals)); % Rand-pos Position depenedent connection probability
    net_rand2_posDepConProb = zeros(size(modelParam.gridxvals)); % Rand-pos Position depenedent connection probability
    distance_counts = zeros(size(modelParam.gridxvals));
    for ithCell = 1:numel(E_indices)
        for jthCell = (ithCell+1):numel(E_indices)
            PFi = PFmat(E_indices(ithCell),:);
            PFj =  PFmat(E_indices(jthCell),:);

            if max(PFi)<analysisParam.minPeakRate || max(PFj)<analysisParam.minPeakRate
                continue
            end

            if analysisParam.singleClust && ...
                    [ sum(network.cluster_mat(:,E_indices(ithCell)))~=1 || ...
                    sum(network.cluster_mat(:,E_indices(jthCell)))~=1 ]
                continue
            end
            if analysisParam.multipleClust && ...
                    [ sum(network.cluster_mat(:,E_indices(ithCell)))<2 || ...
                    sum(network.cluster_mat(:,E_indices(jthCell)))<2 ]
                continue
            end
            %{
            if analysisParam.endClusts && ...
                [ ~[network.cluster_mat(1,E_indices(ithCell))==1 && network.cluster_mat(end,E_indices(ithCell))] ...
                   && ...
                  ~[network.cluster_mat(1,E_indices(jthCell))==1 && network.cluster_mat(end,E_indices(jthCell))] ...
                 ]
                continue
            end
            %}
            nValidPairs = nValidPairs+2; % counts bi-directionally
            [iMax, iInd] = max(PFi);
            [jMax, jInd] = max(PFj);
            if analysisParam.useEV
                iInd = sum(PFi.*[1:numel(PFi)])/sum(PFi);
                jInd = sum(PFj.*[1:numel(PFj)])/sum(PFj);
            end
            pairDist = abs(iInd-jInd)+1;
            if analysisParam.useEV
                pairDist = round(pairDist);
            end
            distance_counts(pairDist) = distance_counts(pairDist)+2; % counts bi-directionally
            pairConns = network.conns(E_indices(ithCell), E_indices(jthCell)) + ...
                network.conns(E_indices(jthCell), E_indices(ithCell));
            net_posDepConProb(pairDist) = net_posDepConProb(pairDist)+ pairConns;
            net_rand_posDepConProb(randi(numel(net_rand_posDepConProb))) = net_rand_posDepConProb(pairDist)+ pairConns;
            net_rand2_posDepConProb(pairDist) = net_rand2_posDepConProb(pairDist) + [rand(1)<0.1];
        end
    end
    net_posDepConProb = net_posDepConProb./nValidPairs;
    net_rand_posDepConProb = net_rand_posDepConProb./nValidPairs;
    net_rand2_posDepConProb = net_rand2_posDepConProb .* [sum(net_posDepConProb)/sum(net_rand2_posDepConProb)];

    op.allPosDepConProb(ithNet,:) = net_posDepConProb;
    op.allPosDepConProb_norm(ithNet,:) = net_posDepConProb./distance_counts*nValidPairs;
    op.allRand_posDepConProb(ithNet,:) = net_rand_posDepConProb;
    op.allRand2_posDepConProb(ithNet,:) = net_rand2_posDepConProb;


    %% Accumulate best events, if there are any
    eventsExist = ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory) && ...
        ~isempty(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue);

    if eventsExist
        rippleEventLengths = resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.endtime-resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.starttime;
        eventRippleIds = find(rippleEventLengths>(modelParam.minEventDur/1000));

        if analysisParam.useWeightedDecode
            rsqrs_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedR2{:});
            rsqrs_shuffle = rsqrs_shuffle(:,analysisParam.ithEnv);
            rsqrs_shuffle = cellfun(@transpose,rsqrs_shuffle,'UniformOutput',false);
            rsqrs_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedR2(:,analysisParam.ithEnv);
            % pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,1); % This is for unweighted pvalue relative to shuffle
            pvals_preplay = nan(size(rsqrs_preplay));
            for ithEvent=1:numel(rsqrs_preplay)
                pvals_preplay(ithEvent) = mean(rsqrs_preplay(ithEvent)<rsqrs_shuffle{ithEvent});
            end
        else
            rsqrs_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_rsquare{:});
            rsqrs_shuffle = rsqrs_shuffle(:,analysisParam.ithEnv);
            rsqrs_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.rsquare(:,analysisParam.ithEnv);

            pvals_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(:,analysisParam.ithEnv); % pvalue relative to shuffle
            %{
            pvals_preplay = nan(size(rsqrs_preplay));
            for ithEvent=1:numel(rsqrs_preplay)
                pvals_preplay(ithEvent) = mean(rsqrs_preplay(ithEvent)<rsqrs_shuffle{ithEvent});
            end
            %}
        end

        % Remove events that have a single non-zero probability time-bin
        improperEventFits = isnan(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.FitpVal(:,analysisParam.ithEnv));
        if analysisParam.removeBadEvents && any(improperEventFits)
            pvals_preplay(improperEventFits)=nan;
            rsqrs_preplay(improperEventFits)=nan;
            rsqrs_shuffle(improperEventFits) = {nan(size(rsqrs_shuffle{1}))};
        end

        maxJumps_preplay = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.maxJump(:,analysisParam.ithEnv);
        maxJumps_shuffle = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_maxJump{:});
        maxJumps_shuffle = maxJumps_shuffle(:,analysisParam.ithEnv);
        op.allEventR2s = [op.allEventR2s; rsqrs_preplay];
        op.allEventPvals = [op.allEventPvals; pvals_preplay];
        op.allEventMaxJumps = [op.allEventMaxJumps; maxJumps_preplay];
        op.allShuffleR2s = [op.allShuffleR2s; rsqrs_shuffle];
        op.allShuffleMaxJumps = [op.allShuffleMaxJumps; maxJumps_shuffle];

        % Extract nTimeBins for each event
        allEventpMats = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{:});
        firstTrajpMats = vertcat(allEventpMats{:,analysisParam.ithEnv});
        firstTrajpMatLengths = arrayfun(@(x) numel(x.timevec), firstTrajpMats);
        op.allEventLengths = [op.allEventLengths; firstTrajpMatLengths];
        op.allShuffleLengths = [op.allShuffleLengths; repelem(firstTrajpMatLengths, numel(rsqrs_shuffle{1}), 1)]; % Note: Shuffles are same length as the corresponding original event

        % Accumulate decode probability across spatial positions
        tmpEvents = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat;
        for ithEvent=1:numel(tmpEvents)
            thisEvent = tmpEvents{ithEvent}{analysisParam.ithEnv}.pMat;
            op.avgDecodePos = op.avgDecodePos + sum(thisEvent, 2)';
            % figure; plot(sum(tmpEvents{ithEvent}{1}.pMat, 2))
            op.allDecodeSums = [op.allDecodeSums, sum(thisEvent, 1)];
        end
    end


    %% Calculate PF score
    if analysisParam.calcScore
        day = 1; epoch = 1; tetrode = 1; tr = 1;
        linfields{day}{epoch}{tetrode}{1}{tr}(:,1) = modelParam.gridxvals*100; % convert to cm
        for ithCell = network.E_indices
            linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5) = PFmat(ithCell,:);
        end
        PFscore = calculate_linfieldsScore(linfields, modelParam, network, trajectory=ithEnv);
        op.PFscores_all = [op.PFscores_all, PFscore];
    end

    % Calculate number of jumps
    %keyboard

    netNJumps = [];
    netNJumpsShuffled = [];
    for ithEvent = 1:numel(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat)
        % Event nJumps
        pMat =  resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{analysisParam.ithEnv}.pMat;
        nEventBins = size(pMat, 2);
        [~, peakBin] = max(pMat);
        decodePeakLocCm = peakBin*modelParam.spatialBin*100;
        eventNJumps = sum(abs(diff(decodePeakLocCm))>analysisParam.jumpThreshCm);
        netNJumps = [netNJumps, eventNJumps];
        % Shuffle nJumps
        for ithShuff = 1:analysisParam.nJumpShuffles
            pMat = pMat(:,randperm(nEventBins));
            [~, peakBin] = max(pMat);
            decodePeakLocCm = peakBin*modelParam.spatialBin*100;

            eventNJumps = sum(abs(diff(decodePeakLocCm))>analysisParam.jumpThreshCm);
            netNJumpsShuffled = [netNJumpsShuffled, eventNJumps];
        end
    end
    op.nJumps = [op.nJumps, netNJumps];
    op.nJumpsShuffled = [op.nJumpsShuffled, netNJumpsShuffled];
    %{
    figure; histogram(netNJumps)
    figure; histogram(netNJumpsShuffled)
    figure; hold on; ecdf(netNJumps); ecdf(netNJumpsShuffled); xline(mean(netNJumps), 'r', LineWidth=2); xline(mean(netNJumpsShuffled))
    [H, P, KSSTAT] = kstest2(netNJumps, netNJumpsShuffled)
    %}

    %% If plotting events...
    if analysisParam.nEventsToPlot>0 && eventsExist
        % PF Peak sequence
        PFmat_E = PFmat(E_indices,:);
        [peakRate,peakRateLocation] = max(PFmat_E, [], 2);
        if analysisParam.useMeanPFDensity % Sort by location of mean PF density
            PFexpectedLocation = sum( PFmat_E./sum(PFmat_E, 2) .*([1:size(PFmat_E, 2)]), 2); % units of space bins
            %PFexpectedLocation(peakRate<analysisParam.minPeakRate)=nan;
            [B,sortedCellIndsbyExpectedLocation] = sort(PFexpectedLocation, 'descend');
            PFpeaksSequence = sortedCellIndsbyExpectedLocation;
        else % Sort by location of PF peak
            %peakRateLocation(peakRate<analysisParam.minPeakRate)=nan;
            [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
            PFpeaksSequence = sortedCellIndsbyPeakRateLocation;
        end

        eventLengths = resultsStruct(ithParam1, ithParam2, ithNet).results.eventLength;
        tmpRanksVec = resultsStruct(ithParam1, ithParam2, ithNet).results.ranksVec;
        eventRelRanks = tmpRanksVec./max(tmpRanksVec);
        eventRelRanks = eventRelRanks(:, eventLengths>=0.05); % temp line, since decoding has different min length
        eventRelRanks(peakRate<analysisParam.minPeakRate,:)=nan; % Exclude non high-rate cells
        [pvals_preplay_sorted, pvals_sorted_inds ]= sort(pvals_preplay);

        % To go from ripple to event:
        % must meet min duration, min cells participating, and
        % participating cells must meet min peak PF rate
        eventRippleIds = utils.getValidEventIDs(PFresultsStruct(ithParam1, ithParam2, ithNet), resultsStruct(ithParam1, ithParam2, ithNet), modelParam);

        % Cluster-wise coloring of raster
        clustMat = network.cluster_mat(:,E_indices)';
        clustColor = clustMat*[1:netParams.clusters]' ./ sum(clustMat, 2);
        % Accumulate best events
        for ithEvent = 1:size(pvals_preplay, 1)
            if any(pvals_preplay(ithEvent) < op.bestEventpVals )
                [~, minInd] = max(op.bestEventpVals);
                op.bestEventpVals(minInd) = pvals_preplay(ithEvent);
                op.bestEventjdVals(minInd) = maxJumps_preplay(ithEvent);
                op.bestEventAbsrVals(minInd) =  sqrt(rsqrs_preplay(ithEvent));
                op.bestEvents_decode{minInd} = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{analysisParam.ithEnv}.pMat;
                op.bestEvents_relRank{minInd} = eventRelRanks(PFpeaksSequence,ithEvent);
                op.bestEvents_relRankColor{minInd} = clustColor(PFpeaksSequence);
                %if isnan(op.bestEventAbsrVals(minInd)); keyboard; end
                %if op.bestEventAbsrVals(minInd)<0.01; keyboard; end

                if analysisParam.loadSpikes
                    if ~isfield(resultsStruct(ithParam1, ithParam2, ithNet).results, 'ripple')
                        warning('Old simulation results structure. No spikes available')
                    end

                    preplayStartInd = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.starttime(eventRippleIds(ithEvent))/modelParam.dt);
                    preplayEndInd = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.endtime(eventRippleIds(ithEvent))/modelParam.dt);
                    eventRaster = netPreplayRaster(:,preplayStartInd:preplayEndInd);
                    eventRaster_E_sort = eventRaster(network.E_indices(PFpeaksSequence),:);
                    op.bestEvents_raster{minInd} = eventRaster_E_sort;

                    runRasterChecks = true;
                    if runRasterChecks
                        % Check that event raster matches results from detect_PBE
                        fracPart = resultsStruct(ithParam1, ithParam2, ithNet).results.eventParticipation(eventRippleIds(ithEvent));
                        rasterActiveCells = find(any( eventRaster, 2 ));
                        rasterActiveCells_E = rasterActiveCells(ismember(rasterActiveCells, network.E_indices));
                        assert(numel(rasterActiveCells_E)/modelParam.n_E==fracPart)

                        % Check that event raster matches results from decode_events()
                        validPlaceCell = utils.getValidPlaceCells(PFresultsStruct(ithParam1, ithParam2, ithNet), modelParam);
                        EcellsSpiked = any(netPreplayRaster(network.E_indices,preplayStartInd:(preplayEndInd-1)), 2);
                        rasterCells = find(any( (EcellsSpiked & validPlaceCell), 2))';
                        decodeCells = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.activecellidx{ithEvent}(:,2)';
                        assert( isequal(rasterCells, decodeCells))
                        disp(['Net ', num2str(ithNet), ', event raster asserts passed'])

                        if 0
                            % Plot current event:
                            MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
                            figure; plotSpikeRaster( [eventRaster(network.E_indices,:); ...
                                eventRaster(network.I_indices,:) ], ...
                                'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
                                'MarkerFormat', MarkerFormat);
                        end
                    end
                end
            end
        end % End loop over net's events

        % Get net's best raster, regardless if it is a best across nets
        if analysisParam.loadSpikes && analysisParam.plotNetsBest && eventsExist
            [pvals_preplay_sorted, pvals_sorted_inds ]= sort(pvals_preplay);
            netBestEventInd = pvals_sorted_inds(1);
            preplayStartInd_best = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.starttime(eventRippleIds(netBestEventInd))/modelParam.dt);
            preplayEndInd_best = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.endtime(eventRippleIds(netBestEventInd))/modelParam.dt);
            eventRaster = netPreplayRaster(:,preplayStartInd_best:preplayEndInd_best);
            eventRaster_E_sort_netBest = eventRaster(network.E_indices(PFpeaksSequence),:);
        end

        op.nEvents = op.nEvents+numel(pvals_preplay);
    end


    %% Net-wise plotting
    switch analysisParam.figSettings
        case 'standard'
            myPlotSettings
        case 'manuscript'
            myPlotSettings(width=1.25, height=1)
        case 'SfNPoster'
            myPlotSettings(width=3, height=2.5, lw=3, afs=24, alw=3)
    end
    % ithNet's structure
    if analysisParam.plotAllNetStruct
        figure; histogram(sum(network.cluster_mat, 1)); xlabel('Clusters (count)'); ylabel('Cells (count)')
        allInputs = vertcat(network.spatialInput{:});
        figure; histogram(allInputs(allInputs~=0)); xlabel('Input strength (pS)'); ylabel('Cells (count)')
        % figure; histogram(sum(network.cluster_mat, 2)); xlabel('n Neurons'); ylabel('Clusters (count)')
    end
    % ithNet's best event
    if analysisParam.plotNetsBest && eventsExist
        tBinSz = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.tBinSz;
        % Plot relative ranks
        %{
        figure; scatter(eventRelRanks(PFpeaksSequence,pvals_sorted(1))', 1:modelParam.n_E, 'k', '|', 'LineWidth', 2 )
        ylabel('Place cell (sorted)'); xlabel('Event relative rank'); title(['Best event of network ', num2str(ithNet)]); ylim([0, modelParam.n_E])
        ylabel('Cell (sorted)'); title ''
        ylim([0, 400]); yticks([0, 200, 400]); %yticklabels({'0', '50', '100'})
        title(['Network ', num2str(ithNet)])
        %}
        if analysisParam.loadSpikes
            %[~,B]=max(eventRaster_E_sort_netBest,[],2);
            %[~,rasterInds]=sort(B);
            %rasterInds = sum(eventRaster_E_sort_netBest, 2)~=0;
            %rasterInds = 1:size(eventRaster_E_sort_netBest,1);
            rasterInds = sum(eventRaster_E_sort_netBest, 2)~=0 & peakRate(PFpeaksSequence)>analysisParam.minPeakRate;

            % Plot raster
            MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
            % figure; plotSpikeRaster( netPreplayRaster, 'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat);
            figure; plotSpikeRaster( eventRaster_E_sort_netBest(rasterInds,:), ...
                'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
                'MarkerFormat', MarkerFormat);
            set(gca,'YDir','normal')
            xlabel('Time (s)'); ylabel('Cell (sorted)');
            title(['Net ', num2str(ithNet), '''s best'])

            % Plot whole raster with random cell sort
            myPlotSettings(width=2.5, height=2)
            MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
            figure; plotSpikeRaster( eventRaster_E_sort_netBest(randperm(size(eventRaster_E_sort_netBest, 1)),:), ...
                'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
                'MarkerFormat', MarkerFormat);
            set(gca,'YDir','normal')
            xlabel('Time (s)'); ylabel('Cell');
            title(['Net ', num2str(ithNet), '''s best'])
        end


        % Plot the event raster sorted by all different trajectories
        if analysisParam.loadSpikes
            %{
            MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
            figure; tiledlayout('flow', TileSpacing='tight', padding='tight');
            set(gca,'YDir','normal')
            xlabel('Time (s)'); ylabel('Cell (sorted)');
            title(['Net ', num2str(ithNet), '''s best'])
            for ithPanelEnv = 1:modelParam.nEnvironments
                nexttile
                
                % Calcualte, PF peak sequence for each environment
                if iscell(PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields)
                    % if modelParam.nEnvironments>1
                    PFmat_thisEnv = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields{ithPanelEnv};
                elseif ismatrix(PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields) && analysisParam.ithEnv==1
                    %else
                    PFmat_thisEnv = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;
                else
                    disp(['PF data error ', num2str(ithParam1), ' ', num2str(ithParam2), ' ', num2str(ithNet)])
                    PFmat_thisEnv = [];
                end
                PFmat_E_thisEnv = PFmat_thisEnv(E_indices,:);
                [peakRate_thisEnv,peakRateLocation_thisEnv] = max(PFmat_E_thisEnv, [], 2);
                if analysisParam.useMeanPFDensity % Sort by location of mean PF density
                    PFexpectedLocation = sum( PFmat_E_thisEnv./sum(PFmat_E_thisEnv, 2) .*([1:size(PFmat_E_thisEnv, 2)]), 2); % units of space bins
                    %PFexpectedLocation(peakRate<analysisParam.minPeakRate)=nan;
                    [B,sortedCellIndsbyExpectedLocation] = sort(PFexpectedLocation, 'descend');
                    PFpeaksSequence_thisEnv = sortedCellIndsbyExpectedLocation;
                else % Sort by location of PF peak
                    %peakRateLocation(peakRate<analysisParam.minPeakRate)=nan;
                    [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation_thisEnv, 'descend');
                    PFpeaksSequence_thisEnv = sortedCellIndsbyPeakRateLocation;
                end
                % TODO: need to re-order eventRaster_E_sort_netBest based on PF peak seq
                rasterInds = sum(eventRaster_E_sort_netBest, 2)~=0 & peakRate_thisEnv(PFpeaksSequence_thisEnv)>analysisParam.minPeakRate;
                plotSpikeRaster( eventRaster_E_sort_netBest(rasterInds,:), ...
                    'TimePerBin', modelParam.dt, 'PlotType', 'scatter', ...
                    'MarkerFormat', MarkerFormat);
            end
            %}
        end

        % Plot decode
        eventPmat = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{netBestEventInd}{analysisParam.ithEnv}.pMat;
        [yt, xt] = size(eventPmat);
        decodeYvals = [1:yt]*(modelParam.spatialBin*100);
        figure; imagesc([1:xt]*(tBinSz), decodeYvals, eventPmat)
        xlabel('Time (ms)'); ylabel('Position (cm)'); title(['Best event of network ', num2str(ithNet)])
        colormap(hot); title ''; colorbar off; caxis(([0, 0.25]))
        yticks([1, 50, 100]); yticklabels({'0', '50', '100'})
        title(['Net ', num2str(ithNet), '''s best'])
        % Add weighted linear correlation to deocde plot
        [X,y] = meshgrid([1:size(eventPmat,2)],[1:size(eventPmat,1)]); w = eventPmat;
        mdl = fitlm(X(:)*tBinSz, y(:)*(modelParam.spatialBin*100),'Weights',w(:));
        %mdl.Coefficients.Estimate(2)
        %mdl.Rsquared.Ordinary
        hold on
        refline(mdl.Coefficients.Estimate(2), mdl.Coefficients.Estimate(1))
        %disp(['Best event of network ', num2str(ithNet)]);
        ylim([1, 101])

        if ~exist('netDecodeColorBarFigHandle', 'var')
            netDecodeColorBarFigHandle = figure;
            title('colorbar', 'fontSize', 19); colorbar; colormap(hot); caxis([0, 0.25])
        end

        % Plot all alternate decodes
        if modelParam.nEnvironments>1
            switch analysisParam.figSettings
                case 'standard'
                    myPlotSettings
                case 'manuscript'
                    myPlotSettings(width=2.5, height=1.75)
                case 'SfNPoster'
                    myPlotSettings(width=3, height=2.5, lw=3, afs=24, alw=3)
            end
            eventPmat = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{netBestEventInd}{analysisParam.ithEnv}.pMat;
            [yt, xt] = size(eventPmat);
            decodeYvals = [1:yt]*(modelParam.spatialBin*100);
            figure; tiledlayout('flow', TileSpacing='tight', padding='tight');
            for ithPanelEnv = 1:modelParam.nEnvironments
                nexttile
                eventPmat = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{netBestEventInd}{ithPanelEnv}.pMat;
                if  mod(ithPanelEnv,2)==0% mod(ithEnvSequence,2)==0 && ~[ithEnvData==ithEnvSequence]
                    % Flip alterante direction sim
                    eventPmat = flipud(eventPmat);
                end
                imagesc([1:xt]*(tBinSz), decodeYvals, eventPmat)
                %xlabel('Time (ms)'); ylabel('Position (cm)'); title(['Best event of network ', num2str(ithNet)])
                colormap(hot); title ''; colorbar off; caxis(([0, 0.25]))
                yticks([1, 50, 100]); yticklabels({'0', '50', '100'})
                % Add weighted linear correlation to deocde plot
                [X,y] = meshgrid([1:size(eventPmat,2)],[1:size(eventPmat,1)]); w = eventPmat;
                mdl = fitlm(X(:)*tBinSz, y(:)*(modelParam.spatialBin*100),'Weights',w(:));
                %mdl.Coefficients.Estimate(2)
                %mdl.Rsquared.Ordinary
                hold on
                refline(mdl.Coefficients.Estimate(2), mdl.Coefficients.Estimate(1))
                ylim([1, 101])

                h = gca; h.XAxis.Visible = 'off';
                if ithPanelEnv==2 || ithPanelEnv==4
                    h = gca; h.YAxis.Visible = 'off';
                end
            end
            disp([
                'Best event shuffle percentiles: ', num2str(100*(1-resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue(netBestEventInd,:)))
                ])
            sgtitle(['Network ', num2str(ithNet)], "fontsize", 10)
        end

        if analysisParam.loadSpikes
            % Plot surrounding population activity
            switch analysisParam.figSettings
                case 'standard'
                    myPlotSettings
                case 'manuscript'
                    myPlotSettings(width=2.5, height=1.75)
                case 'SfNPoster'
                    myPlotSettings(width=3, height=2.5, lw=3, afs=24, alw=3)
            end
            MarkerFormat = struct; MarkerFormat.MarkerSize = 6;
            sortedESpikes = netPreplayRaster(network.E_indices(PFpeaksSequence),:);
            events = [];
            events(:,1) = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.starttime /modelParam.dt);
            events(:,2) = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.endtime   /modelParam.dt);
            t = [0:modelParam.dt:modelParam.t_max_preplay];
            % Plot Raster
            figure;  tiledlayout(2,1, "TileSpacing","tight", "Padding", "compact")
            ax1 = nexttile; hold on
            %{
            if exist('events', 'var')
                for i = 1:size(events, 1)
                    fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, size(sortedESpikes, 1), size(sortedESpikes, 1)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                end
            end
            %}
            plotSpikeRaster( sortedESpikes, ...
                'TimePerBin', modelParam.dt, 'PlotType', 'scatter', 'MarkerFormat', MarkerFormat); %
            ylabel({'Cell', '(sorted)'});
            h = gca; h.XAxis.Visible = 'off';
            set(gca,'YDir','normal')

            % Plot population firing rate
            ax2 = nexttile; hold on
            meanEPopRate = smoothdata(mean(sortedESpikes, 1)/modelParam.dt, 'gaussian', modelParam.PBE_window);
            PBEthresh = max(mean(meanEPopRate)+(modelParam.PBE_zscore*std(meanEPopRate)), modelParam.PBE_min_Hz); % threshold rate for PBE detection
            if exist('events', 'var')
                for i = 1:size(events, 1)
                    fill([t(events(i,:)), fliplr(t(events(i,:)))], [0, 0, ceil(max(meanEPopRate)), ceil(max(meanEPopRate))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
                end
            end
            ylim([0,  ceil(max(meanEPopRate))])
            plot(t, meanEPopRate)
            ylabel({'Pop. rate', '(Hz)'}); xlabel('Time (s)');
            yline(mean(meanEPopRate), '--k') %, {'+0z'})
            %yline(mean(meanEPopRate)+ std(meanEPopRate), '-k', {'+1z'})
            %yline(mean(meanEPopRate)+ 2*std(meanEPopRate), '-k', {'+2z'})
            yline(PBEthresh, '--r', {'Thresh.'}, LabelHorizontalAlignment='right');
            linkaxes([ax1, ax2], 'x')
            xLims = [preplayStartInd_best*modelParam.dt-0.75, preplayEndInd_best*modelParam.dt+0.75];
            xlim(xLims)
            maxInWindow = max( meanEPopRate(:, round(xLims(1)/modelParam.dt): round(xLims(2)/modelParam.dt)));
            ylim([0, maxInWindow*1.1])
        end

    end % End plotting net's best events

    %% Plot each net's PFs for all environments
    if analysisParam.plotNetPFs
        netPFEnvCorrs = nan(modelParam.nEnvironments, modelParam.nEnvironments);
        switch analysisParam.figSettings
            case 'standard'
                myPlotSettings
            case 'manuscript'
                myPlotSettings(width=4, height=4, ttlfsz=1.0)
            case 'SfNPoster'
                myPlotSettings(width=3, height=2.5, lw=3, afs=24, alw=3)
        end
        plotNormRates = true;
        %useMeanPFDensity = false;
        day = 1; epoch = 1; tetrode = 1;

        % Get all PFs
        %[allEnvPFs, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);
        allEnvPFs = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;

        % Get all PF Peak sequences
        PFpeaksSequence = [];
        for ithEnv = 1:modelParam.nEnvironments
            PFmat_E = allEnvPFs{ithEnv}(E_indices,:);
            [peakRate,peakRateLocation] = max(PFmat_E, [], 2);
            if 0%useMeanPFDensity % Sort by location of mean PF density
                PFexpectedLocation = sum( PFmat_E./sum(PFmat_E, 2) .*([1:size(PFmat_E, 2)]), 2); % units of space bins
                %PFexpectedLocation(peakRate<analysisParam.minPeakRate)=nan;
                [B,sortedCellIndsbyExpectedLocation] = sort(PFexpectedLocation, 'descend');
                PFpeaksSequence{ithEnv} = sortedCellIndsbyExpectedLocation;
            else % Sort by location of PF peak
                %peakRateLocation(peakRate<analysisParam.minPeakRate)=nan;
                [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
                PFpeaksSequence{ithEnv} = sortedCellIndsbyPeakRateLocation;
            end
        end

        figure; tiledlayout(modelParam.nEnvironments, modelParam.nEnvironments, 'TileSpacing', 'tight', 'padding', 'tight') %, 'format', 'compact')
        for ithEnv = 1:modelParam.nEnvironments
            linfields{day}{epoch}{tetrode}{1}{ithEnv}(:,1) = modelParam.gridxvals*100; % convert to cm
            posBins = linfields{day}{epoch}{tetrode}{1}{ithEnv}(:,1);
            tmp_PFpeaksSequence = PFpeaksSequence{ithEnv};
            for jthEnv = 1:modelParam.nEnvironments
                PFmat = allEnvPFs{jthEnv};
                PFmat_E = PFmat(network.E_indices,:);
                if  mod(jthEnv,2) ~= mod(ithEnv,2)% mod(ithEnvSequence,2)==0 && ~[ithEnvData==ithEnvSequence]
                    PFmat_E = fliplr(PFmat_E);
                end
                if plotNormRates
                    rateDenom1 = max(PFmat_E(tmp_PFpeaksSequence,:), [], 2);
                    caxmax = 1;
                else
                    rateDenom1 = 1;
                    caxmax = max(PFmat_E, [], 'all');
                end

                % Plot subplot's data
                ithPlot = sub2ind( [modelParam.nEnvironments, modelParam.nEnvironments], ithEnv, jthEnv);
                %subplot(modelParam.nEnvironments, modelParam.nEnvironments, ithPlot);
                nexttile
                imagesc(posBins, 1:modelParam.n_E, PFmat_E(tmp_PFpeaksSequence,:)./rateDenom1 );
                % Put correlation in title
                pixelCorrs = diag(corr(PFmat_E, allEnvPFs{ithEnv}(network.E_indices,:) ));
                % assert(numel(pixelCorrs)==numel(posBins))
                mapCorr = nanmean(pixelCorrs);
                netPFEnvCorrs(ithEnv, jthEnv) = mapCorr;

                if ithEnv~=modelParam.nEnvironments
                    h = gca; h.XAxis.Visible = 'off';
                end
                if jthEnv~=1
                    h = gca; h.YAxis.Visible = 'off';
                end
                %if ithEnv==1 && jthEnv==modelParam.nEnvironments
                %    xlabel('Position (cm)'); ylabel('Cell (column sort)');
                %end
                title(['r=', num2str(mapCorr, '%0.2f')], fontweight='normal')
                if jthEnv==ithEnv
                    nthEnv = ceil(jthEnv/2);
                    if mod(jthEnv,2)==0; trajDir='right';
                    else; trajDir='left'; end
                    title(['Env. ', num2str(nthEnv), ' ', trajDir], fontweight='normal')
                end
            end % jthEnv loop
        end % ithEnv loop
        op.allNetPFEnvCorrs(ithNet,:,:) = netPFEnvCorrs;
    end


    %% Plot each net's PFs for the first environment
    if analysisParam.plotNetPFsEnv1
        switch analysisParam.figSettings
            case 'standard'
                myPlotSettings
            case 'manuscript'
                myPlotSettings(width=1.75, height=1.5)
            case 'SfNPoster'
                %myPlotSettings(width=3, height=2.5, lw=3, afs=24, alw=3)
        end

        % Get all PFs
        %[allEnvPFs, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);
        allEnvPFs = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;

        ithEnv = 1;
        thisNetPFMatE_all = allEnvPFs{ithEnv}(network.E_indices,:);

        row_all_zeros1 = find(all( thisNetPFMatE_all==0, 2)) ;
        row_n_all_zeros1 = find(~all( thisNetPFMatE_all==0, 2)) ;
        [peakRate,peakRateLocation] = max(squeeze(thisNetPFMatE_all(row_n_all_zeros1,:)), [], 2);
        [B,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
        PFpeaksSequence = [row_n_all_zeros1(sortedCellIndsbyPeakRateLocation); row_all_zeros1];
        [peakRate, peakRateLocation_all] = max(thisNetPFMatE_all, [], 2);

        normRates = 1;
        if normRates
            rateDenomAll = max(thisNetPFMatE_all(PFpeaksSequence,:), [], 2);
            caxmax = 1;
        else
            rateDenomAll = 1;
            caxmax = max(thisNetPFMatE_all, [], 'all');
        end

        PFdatatoplot = thisNetPFMatE_all(PFpeaksSequence(rateDenomAll>analysisParam.minPeakRate),:)./rateDenomAll(rateDenomAll>analysisParam.minPeakRate);
        figure; imagesc( (PFdatatoplot(:,5:end)) ); title('All nets'); colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');

        colormap(jet); title ''; colorbar off
        xt = xticks; xticklabels(xt*(modelParam.spatialBin*100)); xlabel('Position (cm)')
        % set(gca,'xtick',[]); xlabel('')
    end

end % ithNet loop
end % op = calcOP();
