function op = calculateEventSpiking(analysisParam, simParam, modelParam, resultsStruct, PFresultsStruct, spikeTimeCell, ithParam1, ithParam2)
% op = grid.calculateEventSpiking(analysisParam, simParam, modelParam, resultsStruct, PFresultsStruct, spikeTimeCell, ithParam1, ithParam2)
%
% Runs the analysis across networks at the specified parameter point and
% returns the results in the output struct op.
%

%% Set up some values
ithTrial = 1; % Assumes a single preplay trial was simulated
assert(size(spikeTimeCell,2)==1)

allParamResults = [resultsStruct(ithParam1, ithParam2, :).results];
maxEventDt = round(max(vertcat(allParamResults.eventLength))./modelParam.dt)+1;

%% Initialize outputs

% Analyses only needed for example param points
if isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) || isequal(analysisParam.paramPoint, "all")
    allCellProbAll = [];
    singleClustProbAll = [];
    multiClustProbAll = [];
end

% Output for spikes/cell/PBE analyses
if analysisParam.runSpikeDistAnalyses
    op.spikesPerPBE = cell(1, modelParam.nNets); % op.spikesPerPBE{ithNet}(ithEvent,ithCell) is a number of spikes in the ithEvent PBE
    op.histFitR2 = [];
end

% Non 'active' based cluster analyses
if analysisParam.runAllClustAnalyses
    spikesPerClustAll = [];
    spikeTimingPerClustAll = [];
    op.nClustMostActive = [];
    op.clustMostActiveDur = [];
end

op.eventClustRMSD = cell(1, modelParam.nNets);
op.eventPropMostActClust = cell(1, modelParam.nNets);

% Output for 'active' cluster analyses
op.clustActiveLengths = [];
op.clustActiveLengths_expanded = [];
op.nClustActive = [];
op.durClustsActive = [];
op.fracMostActive = [];
op.mostToSecondActive = [];
op.coactiveMatrices = zeros(modelParam.nNets, modelParam.clusters, modelParam.clusters);
op.binRMSE = [];
op.corr3Seq = []; % Accumulate booleans for each 3-cluster event if they are correctly ordered
op.corr3SeqRepl = []; % Same as op.corr3Seq but allowing a cluster to be active more than once
op.rsByNClustActive = cell(1, 10);
op.rsByNClustActiveShuffles = cell(1, 10);

op.activationCorrs = [];
op.activationCorrsShuff = [];

% op.adjacentBins = [];
adjacentBinSimilarities = [];
adjacentBinSimilaritiesNorm = [];
adjacentBinSimilaritiesCircShift = [];
adjacentBinSimilaritiesCircShiftNorm = [];
adjacentBinFracClustActive = [];
adjacentBinFracClustActivePreceding = [];


%% Network loop
for ithNet = 1:size(resultsStruct, 3)

    %% Set up before event-loop
    % Recreate network
    E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
    netSeed = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.netSeed;
    netParams = modelParam;
    netParams.(simParam.variedParam(1).name) = simParam.variedParam(1).range(ithParam1);
    netParams.(simParam.variedParam(2).name) = simParam.variedParam(2).range(ithParam2);
    netParams = set_depedent_parameters(netParams);
    network = create_network(netParams, 'seed', netSeed);
    if ~all(network.E_indices==E_indices); error('Incorrect network'); end
    clusterMembershipMatrix = logical(network.cluster_mat(:,network.E_indices));

    % Define variables
    eventRippleIds = utils.getValidEventIDs(PFresultsStruct(ithParam1, ithParam2, ithNet), resultsStruct(ithParam1, ithParam2, ithNet), modelParam);
    nPreplay = size(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pvalue, 1);

    % Convert from spike time cell array to binary spike matrix
    % TODO: convert to utils.spikeTimes2SpikeMat (option "all" or "E" or "I")?
    netPreplaySpikeTimes = squeeze(spikeTimeCell(ithNet,ithTrial,:));
    netPreplayRasterE = zeros(modelParam.n_E, modelParam.t_steps_preplay);
    for ithECell = 1:modelParam.n_E
        ithECellInd = E_indices(ithECell);
        cellSpikeInds = round(netPreplaySpikeTimes{ithECellInd}/modelParam.dt);
        netPreplayRasterE(ithECell,cellSpikeInds) = 1;
    end

    netActivationMats = [];
    netEventRs = [];

    %% Loop over events and get spikes/cell/event
    for ithEvent = 1:nPreplay

        % Get event raster
        preplayStartInd = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.starttime(eventRippleIds(ithEvent))/modelParam.dt);
        preplayEndInd = round(resultsStruct(ithParam1, ithParam2, ithNet).results.ripple.endtime(eventRippleIds(ithEvent))/modelParam.dt);
        eventRasterE = netPreplayRasterE(:,preplayStartInd:preplayEndInd);

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


        %% Example parameter point analyses
        % Calculate mean spikes per dt
        if isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) || isequal(analysisParam.paramPoint, "all")
            nClustParticE = sum(clusterMembershipMatrix, 1);
            eventSingleClustProb = mean(eventRasterE(nClustParticE==1,:), 1);
            eventMultiClustProb = mean(eventRasterE(nClustParticE>1,:), 1);
            eventAllCellProb = mean(eventRasterE, 1);
            if analysisParam.normEventLength
                % Convert event_*ClustProb to a normalized length
                eventSingleClustProb_norm = interp1(1:numel(eventSingleClustProb), eventSingleClustProb, linspace(1, numel(eventSingleClustProb), analysisParam.nEventRateBins));
                eventMultiClustProb_norm = interp1(1:numel(eventMultiClustProb), eventMultiClustProb, linspace(1, numel(eventMultiClustProb), analysisParam.nEventRateBins));
                eventAllCellProb_norm = interp1(1:numel(eventAllCellProb), eventAllCellProb, linspace(1, numel(eventAllCellProb), analysisParam.nEventRateBins));
                % Accumulate to new row
                singleClustProbAll = [singleClustProbAll; eventSingleClustProb_norm];
                multiClustProbAll = [multiClustProbAll; eventMultiClustProb_norm];
                allCellProbAll = [allCellProbAll; eventAllCellProb_norm];
            else
                eventPad = nan(1, maxEventDt-size(eventSingleClustProb, 2) );
                allCellProbAll(end+1,:)       = [eventSingleClustProb,   eventPad];
                singleClustProbAll(end+1,:)   = [eventMultiClustProb,    eventPad];
                multiClustProbAll(end+1,:)    = [eventAllCellProb,       eventPad];
            end
        end


        %% Optional analyses/outputs

        % Spikes/cell/PBE analysis output
        if analysisParam.runSpikeDistAnalyses
            op.spikesPerPBE{ithNet} = [op.spikesPerPBE{ithNet}, sum(eventRasterE, 2)];
        end

        % Cluster-wise spike counts: mean spikes/event for each cluster
        if analysisParam.runAllClustAnalyses
            eventSpikesPerClust = nan(1, modelParam.clusters);
            eventSpikeTimingPerClust = nan(1, modelParam.clusters);
            for ithClust = 1:modelParam.clusters
                clustCells = clusterMembershipMatrix(ithClust,:);
                eventSpikesPerClust(ithClust) = mean(sum(eventRasterE(clustCells,:), 2)); % Mean spikes per PBE out of cells in ithClust
                eventSpikeTimingPerClust(ithClust) = mean(sum(eventRasterE(clustCells,:).*eventTimeVec, 2));
            end
            spikesPerClustAll = [spikesPerClustAll; eventSpikesPerClust];
            spikeTimingPerClustAll = [spikeTimingPerClustAll; eventSpikeTimingPerClust];
        end

        %  Calculate event cluster RMSD
        if analysisParam.runAllClustAnalyses
            expectedAcrossTime = mean(eventClustWiseRateSmoothed, 2);
            eventRmsdAcrossTime = sqrt( mean( (eventClustWiseRateSmoothed-expectedAcrossTime).^2, 2) );
            %figure; plot(event_rmsd_acrossTime); xlabel('Time'); ylabel('RMSD of cluster-wise rates')
            op.eventClustRMSD{ithNet} = [op.eventClustRMSD{ithNet}, mean(eventRmsdAcrossTime)]; % Accumulate
        end

        % Most-active cluster analysis
        % 1) calculate fraction of activity coming from most-active cluster
        % 2) Count # of most active clusters
        if analysisParam.runAllClustAnalyses
            timeNormClustActivity = eventClustWiseRateSmoothed./max(eventClustWiseRateSmoothed, [], 2);
            totalNormClustActivity = sum(timeNormClustActivity, 2);
            propSingleClust = 1./(totalNormClustActivity); % Proportion of activty across time coming from the most active cluster
            wasMostActiveClust = zeros(1, modelParam.clusters);
            clustMostActiveDurEvent = [];
            for ithClust = 1:modelParam.clusters
                A = [timeNormClustActivity(:,ithClust)==1];
                M = reshape(find(diff([0;A;0])~=0), 2, []);
                mostActiveLengths = diff(M);
                [longestMostActiveLength, ~] = max([mostActiveLengths, 0]);
                wasMostActiveClust(ithClust) = longestMostActiveLength>(analysisParam.minDurClustOn/modelParam.dt);
                if longestMostActiveLength>(analysisParam.minDurClustOn/modelParam.dt)
                    validLengths = longestMostActiveLength(longestMostActiveLength>(analysisParam.minDurClustOn/modelParam.dt)) .*modelParam.dt;
                    clustMostActiveDurEvent = [clustMostActiveDurEvent, validLengths];
                end
            end
            nClustMostActive = sum(wasMostActiveClust); % How many clusters were most-active at any time?

            op.clustMostActiveDur = [op.clustMostActiveDur, clustMostActiveDurEvent];
            op.nClustMostActive = [op.nClustMostActive, nClustMostActive]; % Accumulate accross nets
            op.eventPropMostActClust{ithNet} = [op.eventPropMostActClust{ithNet}, mean(propSingleClust)]; % Accumulate

            if any(ismember(analysisParam.exampleEvent, [ithParam1, ithParam2, ithNet, ithEvent], 'rows')) ...
                    && analysisParam.plotExtraFigs
                figure; plot(timeNormClustActivity); xlabel('Time'); ylabel('Cluster rates (norm.)'); ylim([-0.05, 1.05])
                figure; plot(propSingleClust); xlabel('Time'); ylabel('Frac. clust act. from most-act.');
            end
        end


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
        % How long each cluster is active
        durClustsActive = histcounts(mostActiveClust, 1:modelParam.clusters+1).*modelParam.dt;
        % Count of how many different clusters were active
        nClustActive = numel(unique(mostActiveClust));
        % Fraction of activity that came form the most active cluster
        fracMostActive = nanmean(firstToAllActivRatio);

        mostToSecondActive = activityEnvelopes(:,1)./(activityEnvelopes(:,2));
        mostToSecondActive(isinf(mostToSecondActive)) = nan;
        mostToSecondActive = nanmean(mostToSecondActive);

        % Store all active clusters per event, to calculate correlations
        activeClustVec = zeros(1, modelParam.clusters);
        activeClustVec(unique(mostActiveClust)) = 1;
        netActivationMats = [netActivationMats; activeClustVec];

        %% Preplay correlation values split by number of clusters active

        eventR = sqrt(abs(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.weightedR2(ithEvent,analysisParam.ithEnv)));
        netEventRs = [netEventRs; eventR];
        tmp1 = vertcat(resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.shuffle_weightedR2{:});
        tmp2 = tmp1(ithEvent,analysisParam.ithEnv);
        eventShuffleR = sqrt(abs([tmp2{:}]));
        % figure; hold on; xline(eventR); ecdf(eventShuffleR)
        if nClustActive>numel(op.rsByNClustActive)
            disp(['Event with ', num2str(nClustActive), ' active clusters excluded from op.rsByNClustActive'])
        else
            op.rsByNClustActive{nClustActive} = [op.rsByNClustActive{nClustActive}, eventR];
            op.rsByNClustActiveShuffles{nClustActive} = [op.rsByNClustActiveShuffles{nClustActive}, eventShuffleR];
        end


        %% Co-activity matrices
        if ~isempty(mostActiveClust)
            activeClusterSequence = mostActiveClust([true, diff(mostActiveClust')~=0])';
        else
            activeClusterSequence = [];
        end
        for j = 1:numel(activeClusterSequence)-1
            preActiveCluster = activeClusterSequence(j);
            postActiveCluster = activeClusterSequence(j+1);
            op.coactiveMatrices(ithNet, preActiveCluster, postActiveCluster) = ...
                op.coactiveMatrices(ithNet, preActiveCluster, postActiveCluster) + 1/nPreplay;
        end


        %% 3-active cluster sequence

        % If three unique clusters were active, check if they are in the
        % forward or reverse order of their PF biases
        if numel(activeClusterSequence)==3
            eventSeq = activeClusterSequence;
            pfSeq = network.clusterSequence(:,analysisParam.ithEnv);
            [~, eventTrackOrdering] = ismember(eventSeq, pfSeq);
            isFor = issorted(eventTrackOrdering, 'ascend');
            isRev = issorted(eventTrackOrdering, 'descend');
            % corr3SeqRepl allows cluster repeats
            op.corr3SeqRepl = [op.corr3SeqRepl, (isRev|isFor)];
            % corr3Seq doesn't allow cluster repeats
            if numel(unique(activeClusterSequence))==3
                op.corr3Seq = [op.corr3Seq, (isRev|isFor)];
            else
                %disp(eventSeq)
            end
        end

        %% Statistics quantifying bin-wise decode fit

        % Get decode matrix
        eventPMat = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.pMat{ithEvent}{analysisParam.ithEnv}.pMat;
        % Weighted regression on decode matrix
        [X,y] = meshgrid([1:size(eventPMat,2)],[1:size(eventPMat,1)]);
        w = eventPMat;
        mdl = fitlm(X(:),y(:),'Weights',w(:));
        xPred = X(1,:)';
        yPred = predict(mdl, xPred);
        % Deviation from regression
        deviationMatrix = eventPMat.*(y-yPred');
        binRMSE = sqrt( mean(deviationMatrix.^2));
        % Plot
        if any(ismember(analysisParam.exampleEvent, [ithParam1, ithParam2, ithNet, ithEvent], 'rows')) ...
                && analysisParam.plotExtraFigs
            weightedSlope = mdl.Coefficients.Estimate(2);
            weightedR2 = mdl.Rsquared.Ordinary;
            expectedLocation = sum((eventPMat.*y));
            % Plot decode
            figure;
            imagesc(eventPMat)
            hold on
            set(gca,'YDir','normal')
            plot(xPred, yPred, 'r')
            %plot(xPred, expectedLocation, 'k')
            colorbar
            title(['slope=', num2str(weightedSlope), ', R^2=', num2str(weightedR2)], 'FontWeight', 'normal')
            xlabel('Event time bin (10 ms)')
            ylabel('Position bin (2 cm)')
            % Plot deviation of decode probability from linear fit
            figure;
            imagesc(deviationMatrix)
            hold on
            set(gca,'YDir','normal')
            plot(xPred, yPred, 'r')
            %plot(xPred, expectedLocation, 'k')
            colorbar
            title('Decode deviation from fit', 'FontWeight', 'normal')
            xlabel('Event time bin (10 ms)')
            ylabel('Position bin (2 cm)')
            % Plot bin RMSE from linear fit
            figure;
            plot(binRMSE)
            ylim([min([binRMSE, 0]), max(binRMSE)*1.1])
            title('Bin-wise RMSE from fit', 'FontWeight', 'normal')
            xlabel('Event time bin (10 ms)')
            ylabel('RMSE from decode fit')
        end
        % keyboard

        % Determine if a cluster is active and if it is an onset bin
        activeClusterOnsetTime = find(diff(anyClusterIsActive)==1)*modelParam.dt;
        activeClusterOnsetTimeBin = floor(activeClusterOnsetTime*1000/modelParam.tBinSz)+1;
        isActiveOnsetBin = ismember(xPred, activeClusterOnsetTimeBin);
        nBins = numel(xPred);
        nDtPerBin = (modelParam.tBinSz/1000)/modelParam.dt;
        fracBinDtWithActiveCluster = mean(reshape(anyClusterIsActive(1:nDtPerBin*nBins,:), nDtPerBin, nBins), 1);
        fracBinDtWithActiveClusterExpanded = mean(reshape(anyClusterIsActiveExpanded(1:nDtPerBin*nBins,:), nDtPerBin, nBins), 1);

        if any(ismember(analysisParam.exampleEvent, [ithParam1, ithParam2, ithNet, ithEvent], 'rows')) ...
                && analysisParam.plotExtraFigs
            % Plot fracBinDtWithActiveCluster for examle event
            figure;
            plot(fracBinDtWithActiveCluster)
            ylim([0, 1])
            title('Frac. time with active cluster', 'FontWeight', 'normal')
            xlabel('Event time bin (10 ms)')
            ylabel('Frac. decode time bin')
        end

        % Accumulate adjacent-time-bin similarities
        nBins = size(eventPMat, 2);
        eventBinDist = [];
        for ithBin = 1:nBins-1
            p = eventPMat(:,ithBin);
            q = eventPMat(:,ithBin+1);
            distance = utils.emd(p, q);
            eventBinDist = [eventBinDist, distance];
            meanFracActive = mean(fracBinDtWithActiveCluster(ithBin:ithBin+1));
            adjacentBinFracClustActive = [adjacentBinFracClustActive, meanFracActive];
            adjacentBinFracClustActivePreceding = [adjacentBinFracClustActivePreceding, fracBinDtWithActiveCluster(ithBin)];
        end
        % Plot example event
        if any(ismember(analysisParam.exampleEvent, [ithParam1, ithParam2, ithNet, ithEvent], 'rows')) ...
                && analysisParam.plotExtraFigs
            p = eventPMat(:,1);
            q = eventPMat(:,2);
            figure;
            hold on;
            plot(p);
            plot(q);
            xlabel('Pos. bin (2cm)');
            ylabel('Probability');
            title(['Adjacent bins, distance=', num2str(utils.emd(p, q))], 'FontWeight', 'normal')
        end

        % Normalize by all possible pair-wise combinations:
        eventShuffleBinDist = [];
        for ithShuffle = 1:analysisParam.nShuffles
            shuffleInds = randsample(nBins, 2);
            pShuff = eventPMat(:,shuffleInds(1));
            qShuff = eventPMat(:,shuffleInds(2));
            shuffDistance = utils.emd(pShuff, qShuff);
            eventShuffleBinDist = [eventShuffleBinDist, shuffDistance];
        end
        adjacentBinSimilarities = [adjacentBinSimilarities, eventBinDist];
        adjacentBinSimilaritiesNorm = [adjacentBinSimilaritiesNorm, (eventBinDist-nanmean(eventShuffleBinDist))/nanstd(eventShuffleBinDist)];

        % WIP: Similar to above, but circularly shift decodes based on fit
        eventPMatpCircShift = zeros(size(eventPMat));
        baseOffset = 0; % -1*(size(eventPMat, 1)/2 -  mean(yPred)); %size(eventPMat, 1);
        for ithTimeBin  = 1:nBins
            eventPMatpCircShift(:,ithTimeBin) = circshift(eventPMat(:,ithTimeBin), round(baseOffset+yPred(end-ithTimeBin+1)));
        end
        % Plot example shifted event
        if any(ismember(analysisParam.exampleEvent, [ithParam1, ithParam2, ithNet, ithEvent], 'rows')) ...
                && analysisParam.plotExtraFigs
            w = eventPMatpCircShift;
            mdl2 = fitlm(X(:),y(:),'Weights',w(:));
            xPred2 = X(1,:)';
            yPred2 = predict(mdl2, xPred);
            disp(mdl2)
            % Plot
            figure;
            imagesc(eventPMatpCircShift)
            hold on
            set(gca,'YDir','normal')
            plot(xPred2, yPred2, 'r')
            colorbar
            title('Circ. shift to regression', 'FontWeight', 'normal')
            xlabel('Event time bin (10 ms)')
            ylabel('Position bin (2 cm)')

            p = eventPMatpCircShift(:,1);
            q = eventPMatpCircShift(:,2);
            figure;
            hold on;
            plot(p);
            plot(q);
            xlabel('Pos. bin (2cm)');
            ylabel('Probability');
            title(['Adjacent bins, distance=', num2str(utils.emd(p, q))], 'FontWeight', 'normal')
        end

        eventBinDistCircShift = [];
        for ithBin = 1:nBins-1
            pCircShift = eventPMatpCircShift(:,ithBin);
            qCircShift = eventPMatpCircShift(:,ithBin+1);
            distanceCircShift = utils.emd(pCircShift, qCircShift);
            eventBinDistCircShift = [eventBinDistCircShift, distanceCircShift];
        end

        % Shuffle pmat, fit, circle shift, calculate distance between random pair
        if analysisParam.runShuffleFitting
            eventShuffleBinDistCircShift = [];
            for ithShuffle = 1:analysisParam.nShuffles
                % Create shuffle event
                shuffleEventPMat = eventPMat(:,randsample(nBins, nBins));
                % Fit regression
                [X,y] = meshgrid([1:size(shuffleEventPMat,2)],[1:size(shuffleEventPMat,1)]);
                wShuffle = shuffleEventPMat;
                mdl = fitlm(X(:),y(:),'Weights',wShuffle(:));
                xPred = X(1,:)';
                yPred = predict(mdl, xPred);
                shuffleEventPMatpCircShift = zeros(size(shuffleEventPMat));
                baseOffset = 0; %size(eventPMat, 1)/2 - mdl.Coefficients.Estimate(1) - mean(yPred); %size(eventPMat, 1);
                for ithTimeBin  = 1:nBins
                    shuffleEventPMatpCircShift(:,ithTimeBin) = circshift(eventPMat(:,ithTimeBin), round(baseOffset+yPred(end-ithTimeBin+1)));
                end
                shuffleInds = randsample(nBins, 2);
                pShuffCircShift = shuffleEventPMatpCircShift(:,shuffleInds(1));
                qShuffCircShift = shuffleEventPMatpCircShift(:,shuffleInds(2));
                shuffDistanceCircShift = utils.emd(pShuffCircShift, qShuffCircShift);
                eventShuffleBinDistCircShift = [eventShuffleBinDistCircShift, shuffDistanceCircShift];
            end
        else
            eventShuffleBinDistCircShift = nan(size(binRMSE));
        end

        adjacentBinSimilaritiesCircShift = [adjacentBinSimilaritiesCircShift, eventBinDistCircShift];
        adjacentBinSimilaritiesCircShiftNorm = [adjacentBinSimilaritiesCircShiftNorm, (eventBinDistCircShift-nanmean(eventShuffleBinDistCircShift))/nanstd(eventShuffleBinDistCircShift)];
        % keyboard

        % Calculate a "binNormRMSE" where the binRMSE of each bin is
        % normalized by its rmse across shuffles.
        % A given bin might have a high RMSE but low relative to shuffles
        % (e.g. a bin with bimodal probabiliy distribution)
        if analysisParam.runShuffleFitting
            nEventBins = size(eventPMat, 2);
            shuffleBinRMSE = nan(analysisParam.nShuffles, nEventBins);
            [X,y] = meshgrid([1:size(eventPMat,2)],[1:size(eventPMat,1)]);
            for ithShuffle = 1:analysisParam.nShuffles
                % Randomly permute eventPMat columns
                shufflePMat = eventPMat(:, randsample(1:nEventBins, nEventBins));
                % Weighted regression on decode matrix
                w = shufflePMat;
                mdl = fitlm(X(:),y(:),'Weights',w(:));
                xPred = X(1,:)';
                yPred = predict(mdl, xPred);
                % Deviation from regression
                shuffleDeviationMatrix = shufflePMat.*(y-yPred');
                shuffleBinRMSE(ithShuffle,:) = sqrt( mean(shuffleDeviationMatrix.^2));
            end
            binNormRMSE = (binRMSE-mean(shuffleBinRMSE, 1)) ./ std(shuffleBinRMSE, 1);
        else
            binNormRMSE = nan(size(binRMSE));
        end

        % Plot example event figures
        if any(ismember(analysisParam.exampleEvent, [ithParam1, ithParam2, ithNet, ithEvent], 'rows')) ...
                && analysisParam.plotExtraFigs
            figure
            hold on
            plot(binRMSE)
            plot(binNormRMSE)
            title('Example event, RMSE from best fit line')
            xlabel('Event time bin (10 ms wide)')
            ylabel('Probability mass RMSE')
            legend({'Raw', 'Norm to shuffle'}, 'location', 'best')
        end


        %%  Accumulate all event-wise output values
        op.clustActiveLengths = [op.clustActiveLengths, clustActiveLengths];
        op.clustActiveLengths_expanded = [op.clustActiveLengths_expanded, clustActiveLengthsExpanded];
        op.nClustActive = [op.nClustActive, nClustActive];
        op.durClustsActive = [op.durClustsActive, durClustsActive];
        op.fracMostActive = [op.fracMostActive, fracMostActive];
        op.mostToSecondActive = [op.mostToSecondActive, mostToSecondActive];
        % op.binRMSE(:,i): ithParam1, ithParam2, ithNet, ithEvent, IthEnv, ithTimeBin, binRMSE(ithTimeBin)
        op.binRMSE = [op.binRMSE; [xPred, binRMSE', isActiveOnsetBin, fracBinDtWithActiveCluster', fracBinDtWithActiveClusterExpanded', binNormRMSE']];


        %% Verify that the event raster plot matches the expected results saved in resultsStruct
        if analysisParam.runRasterChecks
            % Check that event raster matches results from detect_PBE
            fracPart = resultsStruct(ithParam1, ithParam2, ithNet).results.eventParticipation(eventRippleIds(ithEvent));
            rasterActiveCellsE = find(any( eventRasterE, 2 ));
            %rasterActiveCells_E = rasterActiveCells(ismember(rasterActiveCells, network.E_indices));
            assert(numel(rasterActiveCellsE)/modelParam.n_E==fracPart)
            % Check that event raster matches results from decode_events()
            validPlaceCell = utils.getValidPlaceCells(PFresultsStruct(ithParam1, ithParam2, ithNet), modelParam);
            EcellsSpiked = any(eventRasterE(:,1:end-1), 2);
            rasterCells = find(any( (EcellsSpiked & validPlaceCell), 2))';
            decodeCells = resultsStruct(ithParam1, ithParam2, ithNet).results.replaytrajectory.activecellidx{ithEvent}(:,2)';
            assert( isequal(rasterCells, decodeCells))
            disp(['Net ', num2str(ithNet), ', event raster asserts passed'])
        end


        %% Plot selected example events
        if any(ismember(analysisParam.exampleEvent, [ithParam1, ithParam2, ithNet, ithEvent], 'rows'))
            disp(['Plotting event ', num2str(ithEvent)])
            grid.plotEventSpikingExampleEvent(clusterMembershipMatrix, anyClusterIsActive, eventClustWiseRateSmoothed, eventRasterE, modelParam, analysisParam);
        end


    end % End of loop over events

    %% Network wise cluster co-activation

    %keyboard
    activationCorrs = corr(netActivationMats);
    shuffActivationMats = netActivationMats;
    for i = 1:size(netActivationMats, 2)
        shuffActivationMats(:,i) = circshift(netActivationMats(:,i), randi(size(netActivationMats, 1)));
    end
    shuffActCorrs = corr(shuffActivationMats);
    %figure; imagesc(activationCorrs); colorbar
    %figure; hold on; ecdf(activationCorrs(:)); ecdf(shuffActCorrs(:)); legend({'Actual', 'Shuffle'}, 'location', 'best')
    %[H,P,KSSTAT] = kstest2(activationCorrs(:), shuffActCorrs(:))

    %figure; imagesc(activationCorrs); colorbar


    op.activationCorrs = [op.activationCorrs; activationCorrs(:)];
    op.activationCorrsShuff = [op.activationCorrsShuff; shuffActCorrs(:)];


    %% Post event-loop analyses

    % TODO: turn into function
    % [indexValue, pValue] = utils.calculateSymmetryIndex(X, analysisParam.nShuffles=500);
    % Which index? This one https://math.stackexchange.com/questions/2048817/metric-for-how-symmetric-a-matrix-is ?
    calculateSymmetryIndex = @(x) sqrt(sum( (x-x').^2, 'all' )) / (size(x, 1)^2-size(x, 1));
    xMat = squeeze(op.coactiveMatrices(ithNet,:,:));
    symmetryIndex = calculateSymmetryIndex(xMat);
    shuffleSymmetryIndex = zeros(1, analysisParam.nShuffles);
    for ithShuffle = 1:analysisParam.nShuffles
        xShuffle = xMat(randperm(modelParam.clusters), randperm(modelParam.clusters));
        shuffleSymmetryIndex(ithShuffle) = calculateSymmetryIndex(xShuffle);
    end
    symmetryPVal = mean(symmetryIndex>shuffleSymmetryIndex)';

    op.coactiveSymmetryIndex(ithNet) = symmetryIndex;
    op.coactiveSymmetryPVal(ithNet) = symmetryPVal;

    %% Plot example network figures
    if ithNet==1 && (isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) || isequal(analysisParam.paramPoint, "all")) ...
            if isequal(analysisParam.figSettings, 'manuscript')
            myPlotSettings(width=2.0, height=1.75)
            end
            % Coactivity matrix
            figure;
            imagesc(squeeze(op.coactiveMatrices(ithNet,:,:)));
            cb = colorbar;
            ylabel(cb,'P(c_i followed c_j)') %, 'interpreter', 'latex')
            %title({'Example net, cluster event coactivity', ['Deviation from symmetry=', num2str(symmetryIndex)]})
            title({['Asymmetry index: ', num2str(symmetryIndex)]})
            xlabel('Cluster index, c_i')
            ylabel('Cluster index, c_j')
            % Shuffle histogram
            figure;
            hold on
            xline(symmetryIndex, 'r:', 'lineWidth', 2)
            histogram(shuffleSymmetryIndex);
            title('Example net, cluster event coactivity')
            legend({'Actual', 'Shuffles'}, 'location', 'best')
            ylabel('Shuffles (count)')
            xlabel('Asymmetry index')
    end

    if ithNet==1 && (isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) || isequal(analysisParam.paramPoint, "all")) ...
            && analysisParam.plotExtraFigs
        % Scatter of binRMSE across time bins
        figure;
        scatter(op.binRMSE(:,1), op.binRMSE(:,2))
        title('Example network events')
        xlabel('Event bins')
        ylabel('binRMSE')
        % Histogram of binRMSE
        figure;
        histogram(op.binRMSE(:,2))
        title('Example network, binRMSE')
        xlabel('bin RMSE from fit')
        ylabel('Events (count)')
        % Histogram of binRMSE
        figure;
        histogram(op.binRMSE(op.binRMSE(:,3)==1,2))
        xlabel('binRMSE of isActiveOnsetBin bins')
        ylabel('Events (count)')
    end


end % ithNet loop


% Do events with certain clusters active tend to have higher |r| values?
if 0
    tmp = [];
    for i = 1:modelParam.clusters
        clustIsActive = logical(netActivationMats(:,i));
        meanRActive = mean(netEventRs(clustIsActive));
        meanRNotActive = mean(netEventRs(~clustIsActive));
        meanRAll = mean(netEventRs);
        %tmp = [tmp; (meanRActive / (meanRActive+meanRNotActive)) ];
        tmp = [tmp; (meanRActive / meanRAll) ];
    end
    figure;
    histogram(tmp)
    xlabel('|r| if clust is active (norm)')
    ylabel('Count (clusters)')
    % figure; hold on; ecdf(tmp); xline(1); yline(0.5)
    % figure; scatter(netActivationMats, netEventRs)
end


% figure; hold on; ecdf(op.activationCorrs); ecdf(op.activationCorrsShuff); legend({'Actual', 'Shuffle'}, 'location', 'best')
% [H,P,KSSTAT] = kstest2(op.activationCorrs, op.activationCorrsShuff)
figure; hold on; ecdf(abs(op.activationCorrs)); ecdf(abs(op.activationCorrsShuff));
cdfLegend = legend({'Actual', 'Shuffle'}, 'location', 'best');
cdfLegend.ItemTokenSize = [30*0.5, 18*1];
xlabel('Cluster co-activation |r|')
ylabel({'Cumulative', 'proportion'})
[H,P,KSSTAT] = kstest2(abs(op.activationCorrs), abs(op.activationCorrsShuff));
title(['p=', num2str(P, 2)])

% Example network
figure; imagesc(activationCorrs); colorbar
title('Cluster co-activation correlation')
xlabel('Cluster index'); ylabel('Cluster index')
if 0
    figure; hold on; ecdf(activationCorrs(:)); ecdf(shuffActCorrs(:)); legend({'Actual', 'Shuffle'}, 'location', 'best')
    [H,P,KSSTAT] = kstest2(activationCorrs(:), shuffActCorrs(:))
end


%% Post net-loop analyses

%% Preplay cdfs split by number of clusters active

if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=4, height=2.5)
end

figure;
tiledlayout("flow")
actualCorrVec = [];
actualCorrVecNorm = [];
for iActiveClusts = 1:numel(op.rsByNClustActive)
    if ~isempty(op.rsByNClustActive{iActiveClusts})
        nexttile
        hold on;
        actualRs = op.rsByNClustActive{iActiveClusts};
        shuffleRs = op.rsByNClustActiveShuffles{iActiveClusts};
        medianDiff = median(actualRs) - median(shuffleRs);
        [f1,x1] = ecdf(shuffleRs);
        [f2,x2] = ecdf(actualRs);
        plot(x1, f1, 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
        plot(x2, f2, 'LineWidth', 1, 'color', [0, 0, 0.6])

        actualCorrVec = [actualCorrVec; [repmat(iActiveClusts, numel(actualRs), 1), actualRs']];
        actualCorrVecNorm = [actualCorrVecNorm; [repmat(iActiveClusts, numel(actualRs), 1), [(actualRs-mean(shuffleRs))/std(shuffleRs)]']];

        [~,p_kstest,~] = kstest2(actualRs, shuffleRs);
        disp([ ...
            'ks-test pval for ', ...
            num2str(iActiveClusts), ' active clusters: ', ...
            num2str(p_kstest)...
            ])
        title({[num2str(iActiveClusts), ' active clusters'], ...
            ['median shift ', num2str(medianDiff, 3)]}, ...
            'FontWeight', 'normal', ...
            'FontSize', 10)
        if ismember(iActiveClusts, [1, 4])
            ylabel({'Cumulative', 'Proportion'})
        else
            ylabel('')
        end
        if ismember(iActiveClusts, [4, 5, 6])
            xlabel('|correlation|')
        else
            xlabel('')
        end
    end
end

% Swarm plot of correlations as function of number of Active clusters
if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=2.25, height=1.5)
end
figure;
swarmchart( ...
    actualCorrVecNorm(:,1), ...
    actualCorrVecNorm(:,2), ...
    15, ...
    'blue', ...
    'filled', ...
    'MarkerFaceAlpha',0.25, ...
    'MarkerEdgeAlpha',0.25, ...
    'XJitter','density')
xlabel('Number of active clusters')
ylabel({'|Correlation|', '(z-scored)'})
xticks([unique(actualCorrVecNorm(:,1))])
xlim([min(actualCorrVecNorm(:,1))-0.5, max(actualCorrVecNorm(:,1))+0.5])

mdl = fitlm(actualCorrVecNorm(:, 1), actualCorrVecNorm(:, 2));
hold on;
rf = refline(mdl.Coefficients.Estimate(2), mdl.Coefficients.Estimate(1));
rf.LineStyle = '--';
rf.Color = 'k';
%mdlSlopePVal = mdl.Coefficients.pValue(2);
%mdlr2Val = mdl.Rsquared.Ordinary;

% Spearman rank correlation in title
[rho,pval] = corr(actualCorrVecNorm(:,1), actualCorrVecNorm(:,2), 'Type', 'Spearman');
title(['p=', num2str(pval, 2), ', \rho=', num2str(rho, 2)], 'FontWeight','Normal');


%% Parameter point asymmetry index

if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=2.0, height=1.5)
end
figure;
histogram(op.coactiveSymmetryPVal)
xlabel('Asymmetry index p-value')
ylabel('Networks (count)')


%% Calculate 3-cluster event sequence
if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=2.25, height=1.85)
end

rng(1)

frac3ClustCorr = mean(op.corr3Seq);
nBoots = 500;
replacement = false;
bootSamples = nan(nBoots, numel(op.corr3Seq));
for ithBoot = 1:nBoots
    for ithSeq = 1:numel(op.corr3Seq)
        randSeq = randsample(modelParam.clusters, 3, replacement);
        isFor = issorted(randSeq, 'ascend');
        isRev = issorted(randSeq, 'descend');
        bootSamples(ithBoot, ithSeq) = isFor|isRev;
    end
end
bootMeans = mean(bootSamples, 2);
figure;
title({'Without replacement', ...
    ['Frac 3-seq match: ', num2str(frac3ClustCorr)], ...
    ['percentile: ', num2str(mean(frac3ClustCorr>bootMeans))]}, ...
    'FontWeight', 'Normal')
hold on;
histogram(mean(bootSamples));
xline(frac3ClustCorr, 'r:', 'lineWidth', 2)
xlabel('Matching sequences (frac.)')
ylabel('Shuffles (count)')

frac3ClustCorrRepl = mean(op.corr3SeqRepl);
nBoots = 500;
replacement = true;
bootSamples = nan(nBoots, numel(op.corr3SeqRepl));
for ithBoot = 1:nBoots
    for ithSeq = 1:numel(op.corr3SeqRepl)
        randSeq = randsample(modelParam.clusters, 3, replacement);
        isFor = issorted(randSeq, 'ascend');
        isRev = issorted(randSeq, 'descend');
        bootSamples(ithBoot, ithSeq) = isFor|isRev;
    end
end
bootMeans = mean(bootSamples, 2);
figure;
title({ 'With replacement', ...
    ['Frac 3-seq match: ', num2str(frac3ClustCorrRepl)], ...
    ['percentile: ', num2str(mean(frac3ClustCorrRepl>bootMeans))]}, ...
    'FontWeight', 'Normal')
hold on;
histogram(bootMeans);
xline(frac3ClustCorrRepl, 'r:', 'lineWidth', 2)
xlabel('Matching sequences (frac.)')
ylabel('Shuffles (count)')

disp([num2str(numel(op.corr3SeqRepl) - numel(op.corr3Seq)), ' of ', num2str(numel(op.corr3SeqRepl)), ' 3-cluster events have a cluster repeat'])
disp([num2str(sum(op.corr3SeqRepl)-sum(op.corr3Seq)), ' of those extra events are correctly ordered'])


% Calculate rmsd of spikesPerClust_all
if analysisParam.runAllClustAnalyses
    sortedX = sort(spikesPerClustAll, 2, 'descend');
    eventWiseExpected = mean(spikesPerClustAll, 2);
    eventWiseRMSD = sqrt( mean( (sortedX-eventWiseExpected).^2, 2) );
    op.clustSpikeRMSD = eventWiseRMSD;
end


% Calculate log-linear fit to spikes/cell/PBE distribution
if analysisParam.runSpikeDistAnalyses
    allData = [op.spikesPerPBE{:}];
    [N, edges] = histcounts(allData, 'BinMethod','integers', 'Normalization', 'Probability');
    binCenters = (edges(2:end)+edges(1:end-1)) / 2;
    X = binCenters';
    Y = N';
    Xlog = log(X);
    Ylog = log(Y);
    validInds = ~isinf(Ylog);
    [logLinearFit, gof] = fit(X(validInds),Ylog(validInds),'poly1');
    op.histFitR2 = gof.rsquare;
end


%% Correlation of analyses across adjacent decode bins

if (isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) || isequal(analysisParam.paramPoint, "all")) ...
        && analysisParam.plotExtraFigs

    switch analysisParam.figSettings
        case 'manuscript'
            myPlotSettings(width=3.0, height=2.5)
        otherwise
            myPlotSettings
    end

    figure;
    histogram(adjacentBinSimilarities)
    xlabel('EMD of decoded posteriors')
    ylabel('Adjacent decode bins')

    % Direct comparisons:
    % Fraction of time with active cluster VS decode similarities
    mdl = fitlm(adjacentBinFracClustActive, adjacentBinSimilarities);
    disp(mdl)
    figure;
    plot(mdl)
    legend off
    title(['x0 p=', num2str(mdl.Coefficients.pValue(1)), ', x1 p=', num2str(mdl.Coefficients.pValue(2))])
    xlabel('Active cluster frac. (both bins)')
    ylabel('Adjacent bins EMD')
    % Fraction of time with active cluster VS decode similarities normalized
    x = adjacentBinFracClustActive;
    y = adjacentBinSimilaritiesNorm;
    mdl = fitlm(x, y);
    disp(mdl)
    figure;
    %plot(mdl); legend off
    scatter(x, y, 7.5)
    hline = refline(mdl.Coefficients.Estimate(2),  mdl.Coefficients.Estimate(1));
    hline.Color = 'k';
    hline.LineStyle = "--";
    xlabel('Active cluster frac. (both bins)')
    ylabel('Adjacent bins EMD (norm.)')
    title(...
        { ...
        ['x0 p=', num2str(mdl.Coefficients.pValue(1))], ...
        ['x1 p=', num2str(mdl.Coefficients.pValue(2))] ...
        }, ...
        'FontWeight', 'normal', ...
        'FontSize', 12 ...
        )
    % Circ-shifted comparisons:
    % Fraction of time with active cluster VS decode similarities
    mdl = fitlm(adjacentBinFracClustActive, adjacentBinSimilaritiesCircShift);
    disp(mdl)
    figure;
    plot(mdl)
    legend off
    title(['x0 p=', num2str(mdl.Coefficients.pValue(1)), ', x1 p=', num2str(mdl.Coefficients.pValue(2))])
    xlabel('Active cluster frac. (both bins)')
    ylabel('Adj. bin EMD, circ. adjusted')
    % Fraction of time with active cluster VS decode similarities normalized
    x = adjacentBinFracClustActive;
    y = adjacentBinSimilaritiesCircShiftNorm;
    mdl = fitlm(x, y);
    disp(mdl)
    figure;
    % plot(mdl); legend off
    scatter(x, y, 7.5)
    hline = refline(mdl.Coefficients.Estimate(2),  mdl.Coefficients.Estimate(1));
    hline.Color = 'k';
    hline.LineStyle = "--";
    xlabel('Active cluster frac. (both bins)')
    ylabel('Adj. bin EMD, circ. adjusted (norm.)')
    title(...
        { ...
        ['x0 p=', num2str(mdl.Coefficients.pValue(1))], ...
        ['x1 p=', num2str(mdl.Coefficients.pValue(2))] ...
        }, ...
        'FontWeight', 'normal', ...
        'FontSize', 12 ...
        )
    % Same as above, but x-value is only the precending bins fract cluster active
    %
    mdl = fitlm(adjacentBinFracClustActivePreceding, adjacentBinSimilarities);
    disp(mdl)
    figure;
    plot(mdl)
    legend off
    title(['x0 p=', num2str(mdl.Coefficients.pValue(1)), ', x1 p=', num2str(mdl.Coefficients.pValue(2))])
    xlabel('Active cluster frac. (preceding bin)')
    ylabel('Adjacent bins EMD')
    %
    mdl = fitlm(adjacentBinFracClustActivePreceding, adjacentBinSimilaritiesNorm);
    disp(mdl)
    figure;
    plot(mdl)
    legend off
    title(['x0 p=', num2str(mdl.Coefficients.pValue(1)), ', x1 p=', num2str(mdl.Coefficients.pValue(2))])
    xlabel('Active cluster frac. (preceding bin)')
    ylabel('Adjacent bins EMD (norm.)')
    %
    % Circshift
    %
    mdl = fitlm(adjacentBinFracClustActivePreceding, adjacentBinSimilaritiesCircShift);
    disp(mdl)
    figure;
    plot(mdl)
    legend off
    title(['x0 p=', num2str(mdl.Coefficients.pValue(1)), ', x1 p=', num2str(mdl.Coefficients.pValue(2))])
    xlabel('Active cluster frac. (preceding bin)')
    ylabel('Adj. bin EMD, circ. adjusted')
    %
    mdl = fitlm(adjacentBinFracClustActivePreceding, adjacentBinSimilaritiesCircShiftNorm);
    disp(mdl)
    figure;
    plot(mdl)
    legend off
    title(['x0 p=', num2str(mdl.Coefficients.pValue(1)), ', x1 p=', num2str(mdl.Coefficients.pValue(2))])
    xlabel('Active cluster frac. (preceding bin)')
    ylabel('Adj. bin EMD, circ. adjusted (norm.)')

    keyboard

end


%% Plot parameter point manuscript figures
if (isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) || isequal(analysisParam.paramPoint, "all"))
    if isequal(analysisParam.figSettings, 'manuscript')
        myPlotSettings(width=1.75, height=1.25)
    end
    figure;
    histogram(op.clustActiveLengths_expanded);
    xlabel('Active cluster duration (s)');
    ylabel('Count (active cluster)')
    figure;
    histogram(op.nClustActive);
    xlabel('Active clusters');
    ylabel('Count (events)')
    myPlotSettings
end


%% Plot parameter point extra figures
if (isequal(analysisParam.paramPoint, [ithParam1, ithParam2]) || isequal(analysisParam.paramPoint, "all")) ...
        && analysisParam.plotExtraFigs

    if analysisParam.normEventLength
        xText = 'Time (norm. event length)';
        xLims = [-inf, inf];
    else
        xText = 'Time (dt)';
        xLims = [1, 0.1/modelParam.dt];
    end

    x = singleClustProbAll;
    x_mean = nanmean(x, 1);
    SEM = std(x, [], 1)/sqrt(size(x, 1)); % Standard Error
    ts = tinv([0.025  0.975],size(x, 1)-1); % T-Score
    CI =  x_mean+ ts'*SEM; % Confidence Intervals
    figure;
    hold on;
    plot(x_mean, 'k');
    plot(CI', ':k');
    yline(nanmean(x, 'all'));
    xlabel(xText);
    ylabel('Spike prob. density');
    title('Cells in single cluster')
    xlim(xLims);
    % figure; plot(mean(singleClustProb_all,1)); yline(mean(singleClustProb_all, 'all')); title('Single Cluster')

    x = multiClustProbAll;
    x_mean = nanmean(x, 1);
    SEM = std(x, [], 1)/sqrt(size(x, 1)); % Standard Error
    ts = tinv([0.025  0.975],size(x, 1)-1); % T-Score
    CI =  x_mean+ ts'*SEM; % Confidence Intervals
    figure;
    hold on;
    plot(x_mean, 'k');
    plot(CI', ':k');
    yline(nanmean(x, 'all'));
    xlabel(xText);
    ylabel('Spike prob. density');
    title('Cells in multiple clusters')
    xlim(xLims);
    % figure; plot(mean(multiClustProb_all,1)); yline(mean(multiClustProb_all, 'all')); title('Multi Cluster')

    x = allCellProbAll;
    x_mean = nanmean(x, 1);
    SEM = std(x, [], 1)/sqrt(size(x, 1)); % Standard Error
    ts = tinv([0.025  0.975],size(x, 1)-1); % T-Score
    CI =  x_mean+ ts'*SEM; % Confidence Intervals
    figure;
    hold on;
    plot(x_mean, 'k');
    plot(CI', ':k');
    yline(nanmean(x, 'all'));
    xlabel(xText);
    ylabel('Spike prob. density');
    title('All cells')
    xlim(xLims);
    % figure; plot(mean(singleClustProb_all,1)); yline(mean(singleClustProb_all, 'all')); title('Single Cluster')

    figure;
    plot( nanmean(singleClustProbAll,1) ./ nanmean(allCellProbAll,1) )
    xlabel(xText);
    ylabel('Spike prob (norm. to all cells)');
    title('Cells in single cluster')
    xlim(xLims);

    figure;
    plot( nanmean(multiClustProbAll,1) ./ nanmean(allCellProbAll,1) );
    xlabel(xText);
    ylabel('Spike prob (norm. to all cells)');
    title('Cells in multiple clusters')
    xlim(xLims);

    if analysisParam.runSpikeDistAnalyses
        % log-linear fit to spikes/cell/PBE distribution
        figure;
        hold on
        plot(X, exp(logLinearFit(X)), 'k')
        plot(X, Y, 'c')
        histogram(allData, 'Normalization', 'Probability', 'FaceColor', 'c', 'FaceAlpha', 0.1)
        set(gca,'YScale', 'Log')
        param1Str = [simParam.variedParam(1).name, '=', num2str(simParam.variedParam(1).range(ithParam1))];
        param2Str = [simParam.variedParam(1).name, '=', num2str(simParam.variedParam(2).range(ithParam2))];
        title(['linear log fit: ', param1Str, ', ', param2Str])
        xlabel('Spikes/cell/event'); ylabel('Probability density')

        % Plot PDF fit for example param point
        % Set up data
        allData_shifted = [op.spikesPerPBE{:}]+1;
        [N,edges] = histcounts(allData_shifted, 'BinMethod','integers', 'Normalization', 'Probability');
        binCenters = (edges(2:end)+edges(1:end-1)) / 2;
        X = binCenters';
        Y = N';
        % Fit distributions
        ed = fit(X, Y, 'exp1');
        pd = fit(X, Y, 'power1');
        psd = fit(X, Y, 'x^a*exp(-a)/factorial(round(x))', 'StartPoint', [0.5]);
        % Plot distribution fits
        figure; hold on;
        histogram(allData_shifted, 'Normalization', 'probability' )
        plot(X, feval(ed,X))
        plot(X, feval(pd,X))
        plot(X, feval(psd,X))
        legend({'Data', 'exp1', 'power1', 'poiss'}, 'Location', 'Best')
        ylim([min(N)/10, max(N)*1.5])
        xlim([min(binCenters)-1, max(binCenters)+1])
        set(gca,'YScale', 'Log')
        param1Str = [simParam.variedParam(1).name, '=', num2str(simParam.variedParam(1).range(ithParam1))];
        param2Str = [simParam.variedParam(1).name, '=', num2str(simParam.variedParam(2).range(ithParam2))];
        title(['PDF fit: ' param1Str, ', ', param2Str])
        xlabel('Spikes/cell/event'); ylabel('Probability density')
    end

    if analysisParam.runAllClustAnalyses
        % Cluster spike probability
        figure; plot( mean( sort(spikesPerClustAll, 2, 'descend') ) )
        xlabel('Cluster rank (by event)'); ylabel('Mean spikes (cells in cluster)');
        % ECDFs
        figure; hold on; for ithClust=1:modelParam.clusters; ecdf(spikesPerClustAll(:,ithClust) ); end
        ecdf(spikesPerClustAll(:)); h = gca;
        h.Children(1).Color = 'k';
        h.Children(1).LineWidth = 2;
        xlabel('Mean spikes/event'); ylabel('Cumulative fraction (events)'); title('Each line is a cluster, black is all', 'fontweight', 'normal')

        % Cluster expected spike timing
        figure; plot( mean( sort(spikeTimingPerClustAll, 2, 'ascend') ) )
        xlabel('Cluster rank (by event)'); ylabel('Mean spike timing (cells in cluster)');
        % ECDFs
        figure; hold on; for ithClust=1:modelParam.clusters; ecdf(spikeTimingPerClustAll(:,ithClust) ); end
        ecdf(spikeTimingPerClustAll(:)); h = gca;
        h.Children(1).Color = 'k';
        h.Children(1).LineWidth = 2;
        xlabel('Mean spike timing'); ylabel('Cumulative fraction (events)'); title('Each line is a cluster, black is all', 'fontweight', 'normal')

        % Calculated rmsd of spikesPerClust_all
        figure;
        histogram(eventWiseRMSD);
        xlabel('RMSD cluster-spike count');
        ylabel('Event (count)');

        % Number of clusters meeting the critrion of most-active across events
        figure;
        histogram(op.nClustMostActive);
        xlabel('# clusters active');
        ylabel('Events (count)')

        % Mean (across time bins) of the proportion of activity due to the
        % most-active cluster (might not be normalized as expected)
        figure;
        histogram([op.eventPropMostActClust{:}]);
        xlabel('Prop. activity from most active');
        ylabel('Events (count)')

        % Duration of how long a cluster is the most active (over all periods
        % of all events of all nets
        figure;
        histogram(op.clustMostActiveDur);
        xlabel('Duration most-active (s)');
        ylabel('Periods (count)')
    end

    % "Active" cluster analysis:
    % durClustsActive
    figure;
    histogram(op.durClustsActive(op.durClustsActive~=0));
    xlabel('durClustsActive');
    ylabel('Count (nClusters*nEvents)')
    % durClustsActive, log scale
    figure;
    histogram(op.durClustsActive);
    set(gca,'YScale','Log');
    xlabel('durClustsActive');
    ylabel('Count (nClusters*nEvents)')
    % fracMostActive
    figure;
    histogram(op.fracMostActive);
    xlabel('Frac. activity from most active');
    ylabel('Count (events)')
    % Different between nClustActive and nClustMostActive
    if analysisParam.runAllClustAnalyses
        figure;
        histogram(op.nClustActive-op.nClustMostActive);
        xlabel('nClustActive-nClustMostActive');
        ylabel('Count (events)')
    end

    figure;
    histogram(op.clustActiveLengths);
    xlabel('Active cluster duration (2x, s)');
    ylabel('Count (active cluster)')
end


%% Fraction of total activity from most-active cluster

if isequal(analysisParam.figSettings, 'manuscript')
    myPlotSettings(width=1.75, height=1.25)
end

if 0
    figure;
    hold on
    histogram(op.fracMostActive);
    xline(1/modelParam.clusters, 'r:', 'LineWidth', 2)
    xline(mean(op.fracMostActive), 'k:', 'LineWidth', 2)
    xlabel('Fraction of activity');
    ylabel('Count (events)')
end

figure
hold on
histogram(op.mostToSecondActive)
%xline(nanmean(op.mostToSecondActive), 'k:', 'LineWidth', 2)
xlabel('Activity ratio')
%xlabel('Activity ratio, 1^{st} to 2^{nd}')
ylabel('Count (events)')


end