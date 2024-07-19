%% Analyze the JS single-day W-track decodes
% exptAnalysisPlaceFields.m
%
% Outline of script from exptAnalysisSwrSpiking.m
%
% Runs through all animals

myPlotSettings

%% Set up for analysis

addpath(fullfile("..", "src"))
Config = utils.setConfig;
animalprefixArray = Config.animalPrefixes;
dataDir = Config.exptDataPath;

% Data to analyze
sectionDay = 1; % Always 1, need to fix code below if otherwise
allEpochs = 2:2:16;
allTraj = 1:4;

% Analyze subset of data
%animalprefixArray = {'ER1','JS14'};
allEpochs = 2
% allTraj = 1:2;

% Analysis options
calcPFscores = false
plotExamples = false
calcPeakDiffCorr = false
pfsim.gaussFOLower = [10, 0,  0];   % limits on fit params
pfsim.gaussFOUpper = [30, 100, 10]; % limits on fit params
pfsim.peakTarget   = 15;            % Target peak rate
plotPFs = 0;
analysisParam.minPeakRate = 3; % minimum peak PF rate to be considered a place cell
analysisParam.MSPlots = true;


%% Run analysis, looping over animals

op = nan(10, numel(animalprefixArray), numel(allEpochs), numel(allTraj));
op_allCells = cell(10, numel(animalprefixArray), numel(allEpochs), numel(allTraj));

allPval = [];

all_linfield_mat = []; allLinfieldTraj = 1%[1:4]

for ithAnimal = 1:numel(animalprefixArray)
    animalprefix = animalprefixArray{ithAnimal};
    disp(['Starting animal ', animalprefix])

    animalDir = fullfile(dataDir, [animalprefix, '_direct']);

    [~, hpidx] = decode.getCellIdx(animalDir, animalprefix, sectionDay);

    nCells = size(hpidx, 1);

    %load(sprintf('%s%scellinfo.mat', animalDir, animalprefix))
    load(fullfile(animalDir, [animalprefix, 'linfields01.mat']))
    %load(sprintf('%s%spos01.mat', animalDir, animalprefix))


    %% Epoch loop
    % Plot distribution of number of spikes per replay participating cell

    for ithEpoch = 1:numel(allEpochs)
        sectionEpoch = allEpochs(ithEpoch);

        nTraj = numel( linfields{sectionDay}{sectionEpoch}{hpidx(1, 1)}{hpidx(1, 2)} );

        if plotExamples
            ithCell = 1;
            figure; hold on
            for tr = 1:nTraj
                plot(linfields{sectionDay}{sectionEpoch}{hpidx(ithCell, 1)}{hpidx(ithCell, 2)}{tr}(:,1), linfields{sectionDay}{sectionEpoch}{hpidx(ithCell, 1)}{hpidx(ithCell, 2)}{tr}(:,5) )
            end
            legend({'traj1', 'traj2', 'traj3', 'traj4'}, 'Location', 'Best');
            xlabel('Distance from center arm (m)'); ylabel('Un-norm Lin. Place field')
        end

        pairWiseDists = nan(nTraj, nCells, nCells);
        for ithTraj = 1:numel(allTraj)
            traj = allTraj(ithTraj);

            % Plot all linfields of given animal, epoch, and trajectory
            xvals = linfields{sectionDay}{sectionEpoch}{hpidx(1, 1)}{hpidx(1, 2)}{traj}(:,1);
            yvals = 1:nCells;

            linfield_mat = zeros(nCells, numel(linfields{sectionDay}{sectionEpoch}{hpidx(1, 1)}{hpidx(1, 2)}{traj}(:,1)));
            for ithCell = 1:nCells
                linfield_mat(ithCell,:) =linfields{sectionDay}{sectionEpoch}{hpidx(ithCell, 1)}{hpidx(ithCell, 2)}{traj}(:,5)';
            end

            [peakRate_linfields, peakInd_linfields] = max(linfield_mat, [], 2);
            goodPFs = peakRate_linfields>analysisParam.minPeakRate;
            ePFmat = linfield_mat(goodPFs,:);
            if ismember(traj, allLinfieldTraj)
                all_linfield_mat = [all_linfield_mat; ePFmat(:,1:93)]; % WARNING
            end

            if calcPeakDiffCorr
                pairWiseDists(ithTraj,:,:) = abs(peakInd_linfields-peakInd_linfields');
                pairWiseDists(ithTraj,~goodPFs,~goodPFs) = nan;
            end

            if plotExamples
                [peakRate, peakInd] = max(ePFmat, [], 2);
                [~, cellPeakSort] = sort(peakInd);

                figure; imagesc(xvals, 1:size(ePFmat,1), ePFmat(cellPeakSort, :)./peakRate(cellPeakSort));
                cb = colorbar; ylabel(cb, 'Place fields')
                xlabel('Linear pos. (cm)'); ylabel('Cell ID (sorted)')
            end

            %% Calculate Stats

            % 'MeanRate'
            peakRates = max(ePFmat, [], 2);
            %fprintf('PF mean Peak Rate %0.4f %1.0f %1.0f %1.0f \n', meanPeakRate, ithParam1, ithParam2, ithNet)
            op(1, ithAnimal, ithEpoch, ithTraj) = mean(peakRates);
            op_allCells(1, ithAnimal, ithEpoch, ithTraj) = {peakRates};

            % 'specificity'
            PFmat_exceedesThresh =  double(ePFmat>[0.25*max(ePFmat, [], 2)]); PFmat_exceedesThresh(isnan(ePFmat)) = nan ;
            cellSpecificity  = 1 - nanmean(PFmat_exceedesThresh, 2 );
            %cellSpecificity = 1 - mean( ePFmat> [0.25*max(ePFmat, [], 2)], 2 );


            %fprintf('PF specificity %0.4f %1.0f %1.0f %1.0f \n', mean(cellSpecificity), ithParam1, ithParam2, ithNet)
            op(2, ithAnimal, ithEpoch, ithTraj) = mean(cellSpecificity);
            op_allCells(2, ithAnimal, ithEpoch, ithTraj) = {cellSpecificity};

            % 'information'
            spatialInfo = nanmean( [ePFmat./nanmean(ePFmat, 2)] .* log(( ePFmat+eps )./nanmean(ePFmat, 2) ), 2 );
            %fprintf('PF info. %0.4f %1.0f %1.0f %1.0f \n', nanmean(spatialInfo), ithParam1, ithParam2, ithNet)
            op(3, ithAnimal, ithEpoch, ithTraj) = nanmean(spatialInfo);
            op_allCells(3, ithAnimal, ithEpoch, ithTraj) = {spatialInfo};

            % 'score'
            if calcPFscores && ~isempty(linfield_mat)
                xVals = linfields{sectionDay}{sectionEpoch}{hpidx(1, 1)}{hpidx(1, 2)}{traj}(:,1);
                % Reformat matrix of PFs to struct needed for calculate_linfieldsScore()
                network = struct;
                linfields_for_score = {};
                network.E_indices = 1:nCells;
                network.all_indices = 1:nCells;
                day = 1; epoch = 1; tetrode = 1; tr = 1;
                linfields_for_score{day}{epoch}{tetrode}{1}{tr}(:,1) = xVals; % convert to cm
                for ithCell = network.E_indices
                    linfields_for_score{day}{epoch}{tetrode}{ithCell}{tr}(:,5) = linfield_mat(ithCell,:);
                end

                objScore = calculate_linfieldsScore(linfields_for_score, pfsim, network, trajectory=ithTraj);
                PFscore = -1*objScore; % lower is better for objective score
                fprintf('PF score %0.4f %s %1.0f %1.0f \n', PFscore, animalprefix, sectionEpoch, traj)
                op(4, ithAnimal, ithEpoch, ithTraj) = PFscore;
                op_allCells(4, ithAnimal, ithEpoch, ithTraj) = {PFscore};
            else
                op(4, ithAnimal, ithEpoch, ithTraj) = nan;
                op_allCells(4, ithAnimal, ithEpoch, ithTraj) = {nan};
            end

            %% Calculate place field spatial distribution statistics

            % Select only active, excitatory cells
            %PFinds = [ ismember(1:size(PFmat, 1), E_indices)' & (max(PFmat, [], 2)>=minPeakRate) ];
            %ePFmat = PFmat(PFinds,:);

            nBins = size(ePFmat, 2); % number of spatial bins
            nPCs = size(ePFmat, 1); % number of place cells

            [~, pfPeakSpatialInd] = max(ePFmat, [], 2);
            [peakBinCounts,EDGES] = histcounts(pfPeakSpatialInd, 0.5:1:nBins+0.5);

            % 'linear peaks rmsd'
            op(5, ithAnimal, ithEpoch, ithTraj) = sqrt(mean( (peakBinCounts-(nPCs/nBins)).^2 ) );

            % 'Mean PF rmsd'
            popPF = mean(ePFmat, 1);
            op(6, ithAnimal, ithEpoch, ithTraj) = sqrt(mean( (popPF-(nPCs/nBins)).^2 ) );

            % Frac peaks in center 1/3rd
            firstInd = round(nBins/3*1);
            secondInd = round(nBins/3*2);
            op(7, ithAnimal, ithEpoch, ithTraj) = sum(peakBinCounts(firstInd:secondInd))/nPCs;

            % 'rmsd of the spatial information in each spatial bin'
            % or 'bottom x%-tile of bin spatial info'
            % or 'variance of bin spatial info'
            % or 'mean bin spatial info', this correlates with cell-wise spatial info
            binSpatialInfo = nanmean( [ePFmat./mean(ePFmat, 1)] .* log(( ePFmat+eps )./mean(ePFmat, 1) ), 1 );
            op(8, ithAnimal, ithEpoch, ithTraj) = nanmean(binSpatialInfo);
            %op(9, ithAnimal, ithEpoch, ithTraj) = prctile(binSpatialInfo, 33);
            %op(10, ithAnimal, ithEpoch, ithTraj) = var(binSpatialInfo);
            %op(11, ithAnimal, ithEpoch, ithTraj) = sqrt(mean( (binSpatialInfo-mean(binSpatialInfo)).^2 ) );


        end % Traj loop

        if calcPeakDiffCorr
            xx = reshape(pairWiseDists, nTraj, nCells*nCells);
            assert(isequaln(pairWiseDists(:,1,2), xx(:,2)))
            [rho, pval] = corr(xx', 'rows','complete')
            mdl = fitlm(xx(1,:), xx(3,:));
            figure; plot(mdl);
            title(['Animal ', animalprefix, ', pval: ', num2str(mdl.Coefficients.pValue(2))])
            xlabel('Traj i pairwise peak distance'); ylabel('Traj j pairwise peak distance')
            legend off

            tmp = pval.*(triu(nan(size(pval)), 0)+1);
            allPval = [allPval, tmp(~isnan(tmp))'];
        end

    end % Epoch loop

end % Animal loop


%% Plotting setup

title_list = {'Peak (Hz)', 'Specificity', 'Spatial Info.', 'PF Score (optimal=1)', ...
    'Peak RMSD', 'Mean PF RMSD', 'Frac center 1/3rd', 'binSpatialInfo'};

%{
% Plot a single statistic across epochs
ithState = 3;
X = reshape(op(ithState,:,:,:), prod(size(op, [2, 4])), size(op, 3) );
figure; scatter(1:numel(allEpochs), X)
figure; histogram(X(:))
%}

%% Plot manuscript figures

if analysisParam.MSPlots
    myPlotSettings(width=5.0, height=1.25)
    plotXLims = [0, 30; 0, 1; 0,3];
    binLimits = plotXLims; binLimits(1,1) = analysisParam.minPeakRate;
    nBins = [10, 10, 10];
    figure;
    t = tiledlayout('flow', 'Padding', 'compact'); % 'TileSpacing','compact'
    disp(['Analysis over cells, epoch ', num2str(allEpochs)])
    toPlot = [1:3];
    % toPlot = 1:8;
    for ithStat = toPlot
        X_cells = squeeze(op_allCells(ithStat,:,:,:));
        X = vertcat(X_cells{:});
        ax(ithStat) = nexttile;
        binEdges = linspace(binLimits(ithStat,1), binLimits(ithStat,2), nBins(ithStat)+1);
        histogram(X, binEdges, 'EdgeAlpha', 0.3)
        %{
        if ~isnan(nanmean(X, 'all'))
            xline(nanmean(X, 'all'), 'r')
            yy = ylim;
            text(nanmean(X, 'all'), 0.75*yy(2), num2str(round(nanmean(X, 'all'), 2)), ...
                'FontSize',12, 'Color', 'r')
        end
        %}
        xlabel(title_list{ithStat})
        ylabel('Count')
        xlim(plotXLims(ithStat,:))
        disp(['Range of ithStat ', num2str(min(X)), ' ', num2str(max(X))])
    end
    % linkaxes([ax(:)], 'y')
end

%% Plot all stats, means over PFmaps

myPlotSettings(width=8.5, height=4);
figure
t = tiledlayout('flow', 'Padding', 'compact'); % 'TileSpacing','compact'
title(t, ['Analysis over populations, epoch ', num2str(allEpochs)])
toPlot = [1:3,5:8];
toPlot = 1:8;
for ithStat = toPlot
    X = squeeze(op(ithStat,:,:,:));
    nexttile
    histogram(X, 'EdgeAlpha', 0.3)
    if ~isnan(nanmean(X, 'all'))
        xline(nanmean(X, 'all'), 'r')
        yy = ylim;
        text(nanmean(X, 'all'), 0.75*yy(2), num2str(round(nanmean(X, 'all'), 2)), ...
            'FontSize',12, 'Color', 'r')
    end
    title(title_list{ithStat})
end


%% Plot all stats, without averaging over cells

myPlotSettings(width=6, height=2);
figure;
t = tiledlayout('flow', 'Padding', 'compact'); % 'TileSpacing','compact'
title(t, ['Analysis over cells, epoch ', num2str(allEpochs)])
toPlot = [1:3];
% toPlot = 1:8;
for ithStat = toPlot
    X_cells = squeeze(op_allCells(ithStat,:,:,:));
    X = vertcat(X_cells{:});
    nexttile
    histogram(X, 'EdgeAlpha', 0.3)
    if ~isnan(nanmean(X, 'all'))
        xline(nanmean(X, 'all'), 'r')
        yy = ylim;
        text(nanmean(X, 'all'), 0.75*yy(2), num2str(round(nanmean(X, 'all'), 2)), ...
            'FontSize',12, 'Color', 'r')
    end
    title(title_list{ithStat})
end


%%

myPlotSettings

figure; histogram(allPval)
title(['All correlation p-values, median: ', num2str(median(allPval))])
xlabel('Correlation p-value'); ylabel('All animals, all pairwise trajs')
allPvalMedian = median(allPval)
allPvalMean = mean(allPval)


%%

[peakRate, peakInd] = max(all_linfield_mat, [], 2);
[~, cellPeakSort] = sort(peakInd);

figure; imagesc(xvals, 1:size(all_linfield_mat,1), all_linfield_mat(cellPeakSort, :)./peakRate(cellPeakSort));
colormap(jet)
%cb = colorbar; ylabel(cb, 'Place fields')
xlabel('Linear pos. (cm)'); ylabel('Cell ID (sorted)')


%% Do the different trajectories have different properties?
% op_allCells = cell(10, numel(animalprefixArray), numel(allEpochs), numel(allTraj));

%{
ithStat = 2
for ithTraj1 = 1:4
    %X1_cells = op_allCells(ithStat,:,:,ithTraj1);

    for ithTraj2 = ithTraj1:4
        %X2_cells = op_allCells(ithStat,:,:,ithTraj2);

        [H,P,KSSTAT] = kstest2(vertcat(op_allCells{ithStat,:,:,ithTraj1}), ...
                            vertcat(op_allCells{ithStat,:,:,ithTraj2}) );
        P
    end
end
%}


%% Plot traj-wise histograms

%{
% Plot individual trajectory hist-lines with combined trajectory histogram
figure; hold on; xlabel('Spatial Information (bits)'); ylabel('Place cells (prob.)')
histogram(all_spatialInfo(PCs(:)), EDGES_spatInfo, 'Normalization', normMethod ); xlabel('Spatial Information (bits)'); ylabel('Place cells (prob.)')
for tr = 1:4
	%histogram(spatialInfo(tr, PCs(tr,:)), EDGES_spatInfo, 'Normalization', normMethod, 'FaceColor', color_options(tr,:), 'HandleVisibility','off' ); 
    [counts1, binCenters1] = histcounts(spatialInfo(tr, PCs(tr,:)), EDGES_spatInfo, 'Normalization', normMethod);
    plot(CENTERS_spatInfo, counts1, 'Color', color_options(tr,:));
end
legend({'All traj.', 'traj. 1', 'traj. 2', 'traj. 3', 'traj. 4'})
%}

