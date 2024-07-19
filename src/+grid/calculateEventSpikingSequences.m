function op = calculateEventSpikingSequences(analysisParam, modelParam, paramResultsStruct, paramPFresultsStruct, spikeTimeCell, varargin)
% op = grid.calculateEventSpikingSequences(analysisParam, modelParam, resultsStruct, PFresultsStruct, spikeTimeCell)
%
% Runs the analysis across networks at the specified parameter point and
% returns the results in the output struct op.
%
% Warning: This function assumes the modelParam that has been passed in was
% already updated to match the parameter grid parameter point, e.g.
%    loopModelParam = modelParam;
%    loopModelParam.(simParam.variedParam(1).name) = simParam.variedParam(1).range(param1Ind);
%    loopModelParam.(simParam.variedParam(2).name) = simParam.variedParam(2).range(param2Ind);
%    loopModelParam = set_depedent_parameters(loopModelParam);


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'analysisParam', @isstruct)
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'paramResultsStruct',	@isstruct)
addRequired(inputObj, 'paramPFresultsStruct', @isstruct)
addRequired(inputObj, 'spikeTimeCell', @iscell)
addParameter(inputObj, 'plotFigs', true, @islogical)
addParameter(inputObj, 'matchPFtoShuffle', false, @islogical)
parse(inputObj, analysisParam, modelParam, paramResultsStruct, paramPFresultsStruct, spikeTimeCell, varargin{:});
p = inputObj.Results;


%% Set up some values

% Assumes a single preplay trial was simulated
ithTrial = 1;
assert(size(spikeTimeCell,2)==1)

op = struct;
op.cellClustTuples = []; % accumulates an unknown number of values
op.cellIndepTuples = []; % accumulates values

op.netSpikeRankCorr = []; % Accumulate each event's spike-rank correlations to the place field sequence
op.netSpikeRankCorrShuff = [];
op.netSpikeRankNCells = [];
op.clustSpikeRankCorr = [];
op.clustSpikeRankCorrShuff = [];
op.clustSpikeRankNCells = [];

%% Network loop
for ithNet = 1:modelParam.nNets

    %% Set up before event-loop
    % Recreate network
    netSeed = paramPFresultsStruct(ithNet).results{1}.netSeed;
    network = create_network(modelParam, 'seed', netSeed);
    assert(isequal(network.E_indices, paramPFresultsStruct(ithNet).results{1}.E_indices))
    clusterMembershipMatrix = logical(network.cluster_mat(:,network.E_indices));

    % Define variables
    eventRippleIds = utils.getValidEventIDs(paramPFresultsStruct(ithNet), paramResultsStruct(ithNet), modelParam);
    nPreplay = size(paramResultsStruct(ithNet).results.replaytrajectory.pvalue, 1);

    % Convert from spike time cell array to binary spike matrix
    % TODO: convert to utils.spikeTimes2SpikeMat (option "all" or "E" or "I")?
    netPreplaySpikeTimes = squeeze(spikeTimeCell(ithNet,ithTrial,:));
    netPreplayRasterE = zeros(modelParam.n_E, modelParam.t_steps_preplay);
    for ithECell = 1:modelParam.n_E
        ithECellInd = network.E_indices(ithECell);
        cellSpikeInds = round(netPreplaySpikeTimes{ithECellInd}/modelParam.dt);
        netPreplayRasterE(ithECell,cellSpikeInds) = 1;
    end

    % Get PF data
    PFmat = paramPFresultsStruct(ithNet).results{1}.linfields{analysisParam.ithEnv};
    PFmatE = PFmat(network.E_indices,:);
    [pfPeakRate, pfPeakBin] = max(PFmatE, [], 2);
    nSpatialBins = size(PFmatE, 2);
    switch analysisParam.pfRankMethod
        case 'peak'
            pfSpatialBin = pfPeakBin;
        case 'mean'
            pfSpatialBin = sum((PFmatE./sum(PFmatE, 2)).*[1:nSpatialBins], 2);
        otherwise
            error('Unknown value for analysisParam.pfRankMethod')
    end
    pfSpatialBin(pfPeakRate<analysisParam.minPfRate) = nan;
    if analysisParam.onlySingleClustCells
        pfSpatialBin(sum(clusterMembershipMatrix, 1)~=1) = nan;
    end

    % Convert spatial bin to location on track in cm
    if analysisParam.usePfCmUnits
        pfSpatialLocation = pfSpatialBin*modelParam.spatialBin*100;
    else
        pfSpatialLocation = pfSpatialBin;
    end


    % validPlaceCell = utils.getValidPlaceCells(paramPFresultsStruct(ithNet), modelParam);
    % any([~isnan(pfSpatialBin) & ~validPlaceCell])

    if p.matchPFtoShuffle
        cellidxm = paramResultsStruct(ithNet).results.replaytrajectory.cellidxm;
        shuffledCellidxm = paramResultsStruct(ithNet).results.replaytrajectory.shuffledCellidxm;
        hpnum = size(shuffledCellidxm, 1);
        % Shuffle pfSpatialBin
        unshuffledPfSpatialBin = pfSpatialBin;
        %pfSpatialBin = nan(size(pfSpatialBin));
        for ithIdxmCell = 1:size(cellidxm, 1)
            pfSpatialBin(cellidxm(ithIdxmCell,2)) = unshuffledPfSpatialBin(shuffledCellidxm(ithIdxmCell,2));
        end
    end


    %% Loop over events and get spikes/cell/event

    relRankEventSpikesByCluster = nan(modelParam.n_E, modelParam.clusters, nPreplay);
    relRankEventSpikes = nan(modelParam.n_E, nPreplay);
    relRankBinSpikes = []; % Accumulate

    for ithEvent = 1:nPreplay

        % Get event raster
        preplayStartInd = round(paramResultsStruct(ithNet).results.ripple.starttime(eventRippleIds(ithEvent))/modelParam.dt);
        preplayEndInd = round(paramResultsStruct(ithNet).results.ripple.endtime(eventRippleIds(ithEvent))/modelParam.dt);
        eventRasterE = netPreplayRasterE(:,preplayStartInd:preplayEndInd);

        % Option: invert time for negative sloped events
        decodeSlope = paramResultsStruct(ithNet).results.replaytrajectory.slopes(ithEvent,analysisParam.ithEnv);
        if analysisParam.invertDirectionality && decodeSlope<0
            eventRasterE = fliplr(eventRasterE);
        end

        % Within-cluster relative rank spiking
        for ithClust = 1:modelParam.clusters
            if analysisParam.onlySingleClustCells
                clustECells = (clusterMembershipMatrix(ithClust,:)==1) & (sum(clusterMembershipMatrix, 1)==1);
            else
                clustECells = clusterMembershipMatrix(ithClust,:)==1;
            end
            A = logical(eventRasterE(clustECells,:));
            spikingCellInds = any(A, 2);
            [~, firstSpikeInd] = max(A(spikingCellInds,:),  [], 2);
            normSpikeRank = (firstSpikeInd-min(firstSpikeInd))/(max(firstSpikeInd)-min(firstSpikeInd));
            normSpikeRankFullVec = nan(sum(clustECells), 1);
            normSpikeRankFullVec(spikingCellInds) = normSpikeRank;
            % figure; plotSpikeRaster( logical(eventRasterE(clustECells,:)), 'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
            relRankEventSpikesByCluster(clustECells, ithClust, ithEvent) = normSpikeRankFullVec;

            %keyboard
            % WIP: Within-cluster spike-rank correlation to PF sequence
            clustPfSpatialBin = pfSpatialBin(clustECells);
            goodInds = ~isnan(normSpikeRankFullVec) & ~isnan(clustPfSpatialBin);
            x = normSpikeRankFullVec(goodInds);
            y = clustPfSpatialBin(goodInds);
            if sum(goodInds)>1
                [rho, ~] = corr(x, y, type='Spearman', rows='complete');
                op.clustSpikeRankCorr = [op.clustSpikeRankCorr, abs(rho)];
                [rhoShuff, ~] = corr(x, y(randperm(numel(y))), type='Spearman', rows='complete');
                op.clustSpikeRankCorrShuff = [op.clustSpikeRankCorrShuff, abs(rhoShuff)];
                op.clustSpikeRankNCells = [op.clustSpikeRankNCells, sum(goodInds)];
                % figure; scatter(x, y)
                % figure; scatter(x, y(randperm(numel(y))))
            end

        end
        % keyboard
        % figure; imagesc(relRankEventSpikesByCluster(:,:,ithEvent)); title('Current event relrank by cluster'); colorbar
        % figure; imagesc(relRankEventSpikesByCluster(:,1,ithEvent)); title('Current event relrank for clust1'); colorbar
        % corr(relRankEventSpikesByCluster(:,:,ithEvent)', 'rows','complete')

        % Cluster-independent, within-time-bin relative rank spiking
        eventNDt = (size(eventRasterE, 2)-1);
        nDtPerBin = (modelParam.tBinSz/1000)/modelParam.dt;
        eventNBins = floor(eventNDt*modelParam.dt*1000/modelParam.tBinSz);
        eventBinDtEdges = 1:nDtPerBin:(nDtPerBin*(eventNBins+1));
        if analysisParam.invertDirectionality && decodeSlope<0
            eventDtOverhang = eventNDt-max(eventBinDtEdges)+1;
            eventBinDtEdges = eventBinDtEdges + eventDtOverhang;
        end
        for ithBin = 1:(numel(eventBinDtEdges)-1)
            binSpikes = eventRasterE(:,eventBinDtEdges(ithBin):eventBinDtEdges(ithBin+1));
            %figure; imagesc(binSpikes); title(num2str(ithBin))
            A = logical(binSpikes);
            spikingCellInds = any(A, 2);
            [~, firstSpikeInd] = max(A(spikingCellInds,:),  [], 2);
            normSpikeRank = (firstSpikeInd-min(firstSpikeInd))/(max(firstSpikeInd)-min(firstSpikeInd));
            x = nan(size(spikingCellInds));
            x(spikingCellInds) = normSpikeRank;
            relRankBinSpikes = [relRankBinSpikes, x];
        end

        % Relative rank spiking independent of cluster membership
        A = logical(eventRasterE);
        spikingCellInds = any(A, 2);
        [~, firstSpikeInd] = max(A(spikingCellInds,:),  [], 2);
        normSpikeRank = (firstSpikeInd-min(firstSpikeInd))/(max(firstSpikeInd)-min(firstSpikeInd));
        normSpikeRankFullVec = nan(size(spikingCellInds));
        normSpikeRankFullVec(spikingCellInds) = normSpikeRank;
        relRankEventSpikes(:, ithEvent) = normSpikeRankFullVec;


        % Within-network spike-rank correlation to PF sequence
        goodInds = ~isnan(normSpikeRankFullVec) & ~isnan(pfSpatialBin);
        x = normSpikeRankFullVec(goodInds);
        y = pfSpatialBin(goodInds);
        [rho, ~] = corr(x, y, type='Spearman', rows='complete');
        op.netSpikeRankCorr = [op.netSpikeRankCorr, abs(rho)];
        [rhoShuff, ~] = corr(x, y(randperm(numel(y))), type='Spearman', rows='complete');
        op.netSpikeRankCorrShuff = [op.netSpikeRankCorrShuff, abs(rhoShuff)];
        op.netSpikeRankNCells = [op.netSpikeRankNCells, sum(goodInds)];
        % figure; scatter(x, y)
        % figure; scatter(x, y(randperm(numel(y))))

        %% WIP: Relative-rank analysis within active-event periods
        if 0
            %% Cluster-wise population rates
            % Set up time vector
            nEventDt = size(eventRasterE, 2);
            if analysisParam.normEventLength
                eventTimeVec = 1/nEventDt:1/nEventDt:1;
            else
                eventTimeVec = modelParam.dt:modelParam.dt:(modelParam.dt*nEventDt);
            end

            % Calculate cluster-wise spiking across event time
            eventClustWiseRate = nan(modelParam.clusters, numel(eventTimeVec));
            for ithClust = 1:modelParam.clusters
                clustCells = clusterMembershipMatrix(ithClust,:);
                eventClustWiseRate(ithClust,:) = mean(eventRasterE(clustCells,:), 1) ./ modelParam.dt; % Mean spikes/dt for cells in each cluster
            end
            eventClustWiseRateSmoothed = smoothdata(eventClustWiseRate, 2, 'gaussian', analysisParam.smoothWindow)';

            %% Active cluster analysis
            % Rather than considering the 'most-active' cluster, look for times
            % where any cluster is at least 2x more active than any other

            % activityEnvelopes(:,1) is the activity of the most-active cluster at each time point
            activityEnvelopes = sort(eventClustWiseRateSmoothed, 2, 'descend');
            % How active is the most-active cluster relative to the second most active
            firstToSecondActivRatio = activityEnvelopes(:,1)./activityEnvelopes(:,2);
            firstToAllActivRatio = activityEnvelopes(:,1)./sum(activityEnvelopes(:,1:end), 2);
            anyClusterIsActive = firstToSecondActivRatio>analysisParam.activeClustMuliplyer;

            % Calculate durations of consecutive 1's in anyClusterIsActive
            M = reshape(find(diff([0;anyClusterIsActive;0])~=0),2,[]);
            clustActiveLengths = diff(M)*modelParam.dt; % Durations in (s) of how long any cluster is active

            % Calculate anyClusterIsActiveExpanded as anyClusterIsActive but
            % expanding the periods of active regions until the 2x active
            % cluster becomes no longer the 1x most active cluster
            [~, mostActiveCluster] = max(eventClustWiseRateSmoothed, [], 2);
            anyClusterIsActiveExpanded = anyClusterIsActive;
            % Expand windows to the right
            for ithStep = 2:numel(anyClusterIsActiveExpanded)
                if (anyClusterIsActiveExpanded(ithStep-1)==1) && ...
                        (anyClusterIsActiveExpanded(ithStep)==0)   && ...
                        (mostActiveCluster(ithStep-1)==mostActiveCluster(ithStep))
                    anyClusterIsActiveExpanded(ithStep)=1;
                end
            end
            % Expand windows to the right
            for ithStep = (numel(anyClusterIsActiveExpanded)):-1:2
                if anyClusterIsActiveExpanded(ithStep)==1 && ...
                        anyClusterIsActiveExpanded(ithStep-1)==0   && ...
                        mostActiveCluster(ithStep-1)==mostActiveCluster(ithStep)
                    anyClusterIsActiveExpanded(ithStep-1)=1;
                end
            end
            M = reshape(find(diff([0;anyClusterIsActiveExpanded;0])~=0),2,[]);
            clustActiveLengthsExpanded = diff(M)*modelParam.dt; % Durations in (s) of how long any cluster is active

            % Get index of the clusters that are active during the periods
            % identified by the vector anyClusterIsActive
            [~, mostActiveClust] = max(eventClustWiseRateSmoothed(anyClusterIsActive,:), [], 2);
            % Count of how many different clusters were active
            nClustActive = numel(unique(mostActiveClust));

            % Number of active cluster periods, including possiblity that the
            % same cluster is active more than one
            activeOnsetInds = find(diff([0; anyClusterIsActive])==1);
            nActiveClust = numel(activeOnsetInds);
            % Increment over the periods of active clusters
            for i = 1:numel(nActiveClust)
                tInd = activeOnsetInds(i);
                [~, clustID] = max(eventClustWiseRateSmoothed(tInd,:));
                % Get list of cells that are in the active cluster
                clustECells = clusterMembershipMatrix(clustID,:)==1;
                %figure
                %plotSpikeRaster( logical(eventRasterE(clustECells,:)), ...
                %    'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
                keyboard
                % WIP: calculate relative rank within cluster for each period
            end
        end


        %% Verify that the event raster plot matches the expected results saved in resultsStruct
        if analysisParam.runRasterChecks
            % Check that event raster matches results from detect_PBE
            fracPart = paramResultsStruct(ithNet).results.eventParticipation(eventRippleIds(ithEvent));
            rasterActiveCellsE = find(any( eventRasterE, 2 ));
            %rasterActiveCells_E = rasterActiveCells(ismember(rasterActiveCells, network.E_indices));
            assert(numel(rasterActiveCellsE)/modelParam.n_E==fracPart)
            % Check that event raster matches results from decode_events()
            validPlaceCell = utils.getValidPlaceCells(paramPFresultsStruct(ithNet), modelParam);
            EcellsSpiked = any(eventRasterE(:,1:end-1), 2);
            rasterCells = find(any( (EcellsSpiked & validPlaceCell), 2))';
            decodeCells = paramResultsStruct(ithNet).results.replaytrajectory.activecellidx{ithEvent}(:,2)';
            assert( isequal(rasterCells, decodeCells))
            disp(['Net ', num2str(ithNet), ', event raster asserts passed'])
        end

    end % End of loop over events

    % Post event-loop analyses
    meanRelRankEventSpikesByCluster = nanmean(relRankEventSpikesByCluster, 3);
    meanRelRankEventSpikes = nanmean(relRankEventSpikes, 2);
    meanRelRankBinSpikes = nanmean(relRankBinSpikes, 2);

    % Calculate correlation for just the current network
    if analysisParam.calcCorrForEachNet
        mdl = fitlm(meanRelRankEventSpikes, pfSpatialLocation);
        disp(mdl)
        figure;
        plot(mdl)
        for ithClust = 1:modelParam.clusters
            mdl = fitlm(meanRelRankEventSpikesByCluster(:,ithClust), pfSpatialLocation);
            disp(['Cluster ', num2str(ithClust), ' corr pval ', num2str(mdl.Coefficients.pValue(2))])
        end
        mdl = fitlm(meanRelRankEventSpikesByCluster, pfSpatialLocation);
        disp(mdl)
        % All cluster-pair rel-rank correlations
        relRankCrossClusterCorr = nan(modelParam.clusters);
        for i = 1:modelParam.clusters
            for j = i:modelParam.clusters
                val = corr(meanRelRankEventSpikesByCluster(:,[i,j]), 'rows','complete');
                relRankCrossClusterCorr(i, j) = val(1,2);
                relRankCrossClusterCorr(j, i) = val(1,2);
            end
        end
        figure
        imagesc(relRankCrossClusterCorr, 'AlphaData', ~isnan(relRankCrossClusterCorr))
        colorbar
    end

    % Accumulate all output tuples
    % cellClustTuples: (mean rel cluster rank, PF location)
    % cellIndepTuples: ()
    tmp = [];
    for ithClust = 1:modelParam.clusters
        tmp = [tmp; meanRelRankEventSpikesByCluster(:,ithClust), pfSpatialLocation];
    end
    op.cellClustTuples = [op.cellClustTuples; tmp(~any(isnan(tmp), 2), :)];
    op.cellIndepTuples = [op.cellIndepTuples; meanRelRankEventSpikes, pfSpatialLocation, meanRelRankBinSpikes];

    % Plot example figures
    if (ithNet==1 || analysisParam.plotExtraFigs) && p.plotFigs
        % Mean rel rank for all clusters
        plotAllRelRankByCluster(meanRelRankEventSpikesByCluster, ithNet);
        % Mean+/-SEM for example cluster
        plotExampleClusterRelRank(relRankEventSpikesByCluster, ithNet);
        % Plot mean-rel rank independent of clusters
        plotAllCellRelRank(relRankEventSpikes, ithNet);
    end

end % ithNet loop


%% Analyses across all networks

% Event spike rank correlations
if p.plotFigs
    if isequal(analysisParam.figSettings, 'manuscript')
        myPlotSettings(width=2.5, height=1.5)
    end
    boolDict = dictionary(0, 'False', 1, 'True');

    % Within-network PBE spiking vs PF location rank correlations
    [~,p_kstest,~] = kstest2(op.netSpikeRankCorr, op.netSpikeRankCorrShuff);
    figure;
    hold on;
    [f1,x1] = ecdf(op.netSpikeRankCorrShuff);
    [f2,x2] = ecdf(op.netSpikeRankCorr);
    plot(x1, f1, 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
    plot(x2, f2, 'LineWidth', 1, 'color', [0, 0, 0.6])
    %legend({'Shuffles', 'Preplays'}, 'Location', 'Best', 'Box', 'off')
    title({ ...
        'Within-network PBE spiking vs PF' ...
        ['Single-cluster cells='+boolDict(analysisParam.onlySingleClustCells)] ...
        ['p=', num2str(p_kstest)]}, ...
        fontsize=10, FontWeight='Normal')
    ylabel({'Cumulative', 'Proportion'})
    xlabel('|rank correlation|')

    % Within-cluster PBE spiking vs PF location rank correlations
    [~,p_kstest,~] = kstest2(op.clustSpikeRankCorr, op.clustSpikeRankCorrShuff);
    figure;
    hold on;
    [f1,x1] = ecdf(op.clustSpikeRankCorrShuff);
    [f2,x2] = ecdf(op.clustSpikeRankCorr);
    plot(x1, f1, 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
    plot(x2, f2, 'LineWidth', 1, 'color', [0, 0, 0.6])
    %legend({'Shuffles', 'Preplays'}, 'Location', 'Best', 'Box', 'off')
    title({ ...
        'Within-cluster PBE spiking vs PF' ...
        ['Single-cluster cells='+boolDict(analysisParam.onlySingleClustCells)] ...
        ['p=', num2str(p_kstest)]}, ...
        fontsize=10, FontWeight='Normal')
    ylabel({'Cumulative', 'Proportion'})
    xlabel('|rank correlation|')

    % Scatter of nParticipatingCells with |corr|
    %figure; scatter(op.netSpikeRankCorr, op.netSpikeRankNCells)
    %figure; scatter(op.clustSpikeRankCorr, op.clustSpikeRankNCells)
    %figure; plot(fitlm(op.netSpikeRankCorr, op.netSpikeRankNCells))
    %figure; plot(fitlm(op.clustSpikeRankCorr, op.clustSpikeRankNCells))

    % Removing few-cell events removes the |corr| cdf jumps
    %{
    figure; hold on;
    ecdf(op.clustSpikeRankCorr(op.clustSpikeRankNCells>10))
    ecdf(op.clustSpikeRankCorrShuff(op.clustSpikeRankNCells>10))
    %}

end
if p.plotFigs
    % Cluster-dependent
    [~, op.clustMdl] = plotScatterCorr(op.cellClustTuples(:,2), op.cellClustTuples(:,1), analysisParam);
    xlabel('PF location (cm)')
    ylabel({'Within-cluster', 'Mean relative rank'})
    %ylabel({'Mean rel. within-cluster', 'intra-event rank'})
    disp(['clustMdl: slope=', num2str(op.clustMdl.Coefficients.Estimate(2)), ' p=', num2str(op.clustMdl.Coefficients.pValue(2))])

    % Cluster-independent
    [~, op.netMdl] = plotScatterCorr(op.cellIndepTuples(:,2), op.cellIndepTuples(:,1), analysisParam);
    xlabel('PF location (cm)')
    ylabel({'Within-network', 'Mean relative rank'})
    %ylabel({'Mean rel. cluster-indep.', 'intra-event rank'})
    disp(['netMdl: slope=', num2str(op.netMdl.Coefficients.Estimate(2)), ' p=', num2str(op.netMdl.Coefficients.pValue(2))])

    [~, op.binMdl] = plotScatterCorr(op.cellIndepTuples(:,2), op.cellIndepTuples(:,3), analysisParam);
    xlabel('PF location (cm)')
    ylabel({'Within-bin', 'Mean relative rank'})
    %ylabel({'Mean rel. cluster-indep.', 'intra-bin rank'})
    disp(['binMdl: slope=', num2str(op.binMdl.Coefficients.Estimate(2)), ' p=', num2str(op.binMdl.Coefficients.pValue(2))])
else
    % Cluster-dependent
    op.clustMdl = fitlm(op.cellClustTuples(:,2), op.cellClustTuples(:,1));
    % Cluster-independent
    op.netMdl = fitlm(op.cellIndepTuples(:,2), op.cellIndepTuples(:,1));
    % Within bin sequences
    op.binMdl = fitlm(op.cellIndepTuples(:,2), op.cellIndepTuples(:,3));
end

end % End function

%% Helper functions

function [figHandle, mdl] = plotScatterCorr(x, y, analysisParam)
mdl = fitlm(x, y);
%disp(mdl)
figHandle = figure;
scatter(x, y, 7.5)
hline = refline(mdl.Coefficients.Estimate(2),  mdl.Coefficients.Estimate(1));
hline.Color = 'k';
hline.LineStyle = "--";
% h = plot(mdl); h(1).Marker = 'o'; h(1).MarkerSize = 5; h(1).LineWidth = 0.25;; legend off
boolDict = dictionary(0, 'False', 1, 'True');
title( ...
    { ...
    ['Invert-reverse events='+boolDict(analysisParam.invertDirectionality)] ...
    ['Single-cluster cells='+boolDict(analysisParam.onlySingleClustCells)] ...
    ['Slope pval: ', num2str(mdl.Coefficients.pValue(2))] ...
    }, ...
    'FontWeight', 'normal', ...
    'FontSize', 12 ...
    )
end

function plotAllRelRankByCluster(meanRelRankEventSpikesByCluster, ithNet)
figure;
imagesc(meanRelRankEventSpikesByCluster, 'AlphaData', ~isnan(meanRelRankEventSpikesByCluster))
colorbar;
xlabel('Cluster index')
ylabel('Cell index')
title( ...
    {['Net ', num2str(ithNet),], 'mean rel. within-cluster event rank'}, ...
    'FontWeight', 'normal', ...
    'FontSize', 12 ...
    )
end

function plotExampleClusterRelRank(relRankEventSpikesByCluster, ithNet)
ithCluster = 10;
stdRelRankEventSpikesByCluster = nanstd(relRankEventSpikesByCluster, [], 3);
nCellEvents = sum(~isnan(relRankEventSpikesByCluster), 3);
meanRelRankEventSpikesByCluster = nanmean(relRankEventSpikesByCluster, 3);
y = meanRelRankEventSpikesByCluster(:,ithCluster);
ye = stdRelRankEventSpikesByCluster(:,ithCluster)./sqrt(nCellEvents(:,ithCluster));
x = 1:sum(~isnan(y));
y = y(~isnan(y));
ye = ye(~isnan(ye));
%figure; errorbar(x, y, ye, 'LineStyle','none','Marker','.', 'MarkerSize',10)
%
figure;
[~, sortInds] = sort(y);
errorbar(x, y(sortInds), ye(sortInds), 'LineStyle','none','Marker','.', 'MarkerSize',10)
yline(0.5, "LineStyle", ":")
xlabel('Cell index (sorted)');
ylabel({'Mean rel. within-cluster', 'event rank (\pm SEM)'});
title(['Net ', num2str(ithNet), ', Cluster ', num2str(ithCluster)], ...
    'FontWeight', 'normal', ...
    'FontSize', 12 ...
    )
end

function plotAllCellRelRank(relRankEventSpikes, ithNet)
stdRelRankEventSpikes = nanstd(relRankEventSpikes, [], 2);
nCellEvents = sum(~isnan(relRankEventSpikes), 2);
meanRelRankEventSpikes = nanmean(relRankEventSpikes, 2);
y = meanRelRankEventSpikes;
ye = stdRelRankEventSpikes./sqrt(nCellEvents);
x = 1:sum(~isnan(y));
y = y(~isnan(y));
ye = ye(~isnan(ye));
%figure; errorbar(x, y, ye, 'LineStyle','none','Marker','.', 'MarkerSize',10)
[~, sortInds] = sort(y);
figure;
errorbar(x, y(sortInds), ye(sortInds), 'LineStyle','none','Marker','.', 'MarkerSize',10, 'CapSize', 0)
yline(0.5, "LineStyle", ":")
xlabel('Cell index (sorted)');
ylabel({'Mean rel. cluster-indep.', 'event rank (\pm SEM)'});
title(['Net ', num2str(ithNet)], ...
    'FontWeight', 'normal', ...
    'FontSize', 12 ...
    )
xlim([0, numel(y)])
end
