%% Calculates and plots several statistics of simulated Place fields
% gridAnalysisPlaceFieldStats.m
%
% Meak peak rate, KS-stat, PF specificity, PF spatial info., etc.


%% Choose which simulation results to analyze

decodeFileID = "mainPaperGrid"; % Select which decoding file result ID to analyze


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


%% Analysis parameters and options
analysisParam.paramPoint = [2, 3]; % Select which example param point to plot
analysisParam.ithEnv = 1; % which environments decoding and place fields to plot/analyze
analysisParam.minPeakRate = 3; % minimum peak PF rate to be considered a place cell

% Plotting options
analysisParam.MSPlots = true;
analysisParam.plotExtraFigs = false; % Only plot manuscript figures if false
analysisParam.plotExampleParamPoint = true;
analysisParam.plotPFs = false; % if true and analysisParam.calcPFscores=true, plot place fields of every network
analysisParam.calcPFscores = false;

analysisParam.useParfor = false; % Useful for calcPFscores, which is slow

analysisParam.excludeFullConn = true;
analysisParam.fullConnBoundary  = @(mnc) (mnc.^2)/modelParam.conn_prob;

disp(['Decoding file: ', decodeName])
disp(analysisParam)
disp({simParam.variedParam(:).name})
disp(['Example Param vals.: ', ...
    num2str(simParam.variedParam(1).range( analysisParam.paramPoint(1))), ...
    ' ', num2str(simParam.variedParam(2).range( analysisParam.paramPoint(2))) ])


%% Initialize the output structure op with desired fields, using the format:
%   op.meanPFPeaks(ithParam1, ithParam2, ithNet)

PFStatFields = {'meanPFPeaks', 'meanIRate', 'meanCellSpecificity', 'meanSpatialInfo'};
PFDistFields = {'linearPeakRMSD', 'meanPFRMSD', 'fracPeakCenter3rd', 'spatialBinInfoMean', ...
    'spatialBinInfo33pct', 'spatialBinInfoVar', 'spatialBinInfoRMSD', 'peakDistKLD'};

opFields = [PFStatFields, PFDistFields, {'PFscore', 'exampleCells'}];
tempCell = cell(length(opFields),1);
tempCell(:)= {nan(numel(simParam.variedParam(1).range), numel(simParam.variedParam(2).range),simParam.nNets)};
op = cell2struct(tempCell,opFields);


%% Loop over all parameter sets and all networks

rng('default')
tic
op = loopOPFields(simParam, analysisParam, opFields, PFresultsStruct, modelParam);
runTime = toc;
disp([ 'Analysis runTime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS') ])


%% Plotting setup
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

% Matrix of values to exclude if analysisParam.excludeFullConn
if analysisParam.excludeFullConn && isequal(simParam.variedParam(1).name, 'mnc') && isequal(simParam.variedParam(2).name, 'clusters')
    mncVec = simParam.variedParam(1).range;
    nClustVec = simParam.variedParam(2).range;
    excludedBadConnProb = [(nClustVec' > analysisParam.fullConnBoundary(mncVec))]';
else
    excludedBadConnProb = zeros(numel(simParam.variedParam(1).range), numel(simParam.variedParam(2).range));
end

%% Manuscript figures
if analysisParam.MSPlots

    % PF statistics across parameter grid
    myPlotSettings(width=6.0, height=2.5, ttlfsz=1.0)
    %plotFields = {'meanPFPeaks',  'meanCellSpecificity', 'meanSpatialInfo',     'meanIRate',        'linearPeakRMSD',  'fracPeakCenter3rd'};
    %plotTitles = {'Peak (Hz)',    'Specificity',         'Spatial info.',       'Mean I rate (Hz)', 'Peak dist. RMSE', 'Frac. peaks center'};
    plotFields = {'meanPFPeaks',  'meanCellSpecificity', 'meanSpatialInfo',     'meanIRate',        'peakDistKLD',  'fracPeakCenter3rd'};
    plotTitles = {'Peak (Hz)',    'Specificity',         'Spatial info.',       'Mean I rate (Hz)', 'Peak dist. KL-D', 'Frac. peaks center'};
    if all(cellfun(@(x) any(strcmp(x, opFields)), plotFields))
        figure
        for ithField = 1:numel(plotFields)
            currentField = plotFields{ithField};
            ax = subplot(2,3,ithField);
            imagesc(xParamvec, yParamvec, mean([op.(currentField)], 3)', 'AlphaData', ~isnan(mean([op.(currentField)],3)')&~excludedBadConnProb')
            set(gca,'YDir','normal')
            cb = colorbar(); cb.Label.String = '';
            %xlabel(xName,'Interpreter','none')
            %ylabel(yName,'Interpreter','none')
            title(plotTitles{ithField}, 'FontWeight','Normal')
            if any(strcmp(currentField, {'linearPeakRMSD', 'meanPFRMSD', 'peakDistKLD'} ))
                colormap(ax, flipud(parula))
            end
        end
    else
        warning('Missing desired field in op for plot')
    end

    % PF stat distributions over cells at particular parameter point
    myPlotSettings(width=5.0, height=1.25)
    plotFields = {'Peak (Hz)', 'Specificity', 'Spatial Info.'};
    plotXLims = [0, 30; 0, 1; 0,3];
    %plotXLims = [0, 60; 0, 1; 0,3];
    binLimits = plotXLims; binLimits(1,1) = analysisParam.minPeakRate;
    nBins = [10, 10, 10];
    figure;
    t = tiledlayout('flow', 'Padding', 'compact'); % 'TileSpacing','compact'
    for ithStat = 1:numel(plotFields)
        ithStatData = vertcat( op.exampleCells{ithStat,:,:,:});
        nexttile
        binEdges = linspace(binLimits(ithStat,1), binLimits(ithStat,2), nBins(ithStat)+1);
        histogram(ithStatData, binEdges, 'EdgeAlpha', 0.3)
        %{
            if ~isnan(nanmean(ithStatData, 'all'))
                xline(nanmean(ithStatData, 'all'), 'r')
                yy = ylim;
                text(nanmean(ithStatData, 'all'), 0.75*yy(2), num2str(round(nanmean(X, 'all'), 2)), ...
                    'FontSize',12, 'Color', 'r')
            end
        %}
        xlabel(plotFields{ithStat})
        ylabel('Count')
        xlim(plotXLims(ithStat,:))
        disp(['Range of ithStat ', num2str(min(ithStatData)), ' ', num2str(max(ithStatData))])
    end
end


%% Plot PF statistics

if analysisParam.plotExtraFigs || ~analysisParam.MSPlots
    myPlotSettings(width=4.25, height=3.5)

    plotFields = {'meanPFPeaks', 'meanIRate', 'meanCellSpecificity', 'meanSpatialInfo'};
    if all(cellfun(@(x) any(strcmp(x, opFields)), plotFields))
        figure
        for ithField = 1:numel(plotFields)
            currentField = plotFields{ithField};
            ax = subplot(2,2,ithField);
            imagesc(xParamvec, yParamvec, mean([op.(currentField)], 3)', 'AlphaData', ~isnan(mean([op.(currentField)],3)')&~excludedBadConnProb')
            set(gca,'YDir','normal')
            cb = colorbar(); cb.Label.String = '';
            xlabel(xName,'Interpreter','none')
            ylabel(yName,'Interpreter','none')
            title(currentField)
        end
    else
        warning('Missing desired field in op for plot')
    end
end

%% Plot place field spatial distribution

if analysisParam.plotExtraFigs || ~analysisParam.MSPlots
    plotFields = {'linearPeakRMSD', 'meanPFRMSD', 'fracPeakCenter3rd', 'spatialBinInfoMean'};
    if all(cellfun(@(x) any(strcmp(x, opFields)), plotFields))
        figure
        for ithField = 1:numel(plotFields)
            currentField = plotFields{ithField};
            ax = subplot(2,2,ithField);
            imagesc(xParamvec, yParamvec, mean([op.(currentField)], 3)', 'AlphaData', ~isnan(mean([op.(currentField)], 3)' )&~excludedBadConnProb')
            set(gca,'YDir','normal')
            cb = colorbar(); cb.Label.String = '';
            xlabel(xName,'Interpreter','none')
            ylabel(yName,'Interpreter','none')
            title(currentField)
            if any(strcmp(currentField, {'linearPeakRMSD', 'meanPFRMSD'} ))
                colormap(ax, flipud(parula))
            end
        end
    else
        warning('Missing desired field in op for plot')
    end
end

%% Plot additional Bin-wise spatial info stats

if analysisParam.plotExtraFigs || ~analysisParam.MSPlots
    plotFields = {'spatialBinInfoMean', 'spatialBinInfo33pct', 'spatialBinInfoVar', 'spatialBinInfoRMSD'};
    if all(cellfun(@(x) any(strcmp(x, opFields)), plotFields))
        figure
        for ithField = 1:numel(plotFields)
            currentField = plotFields{ithField};
            ax = subplot(2,2,ithField);
            imagesc(xParamvec, yParamvec, mean([op.(currentField)], 3)', 'AlphaData', ~isnan(mean([op.(currentField)], 3)' )&~excludedBadConnProb')
            set(gca,'YDir','normal')
            cb = colorbar(); cb.Label.String = '';
            xlabel(xName,'Interpreter','none')
            ylabel(yName,'Interpreter','none')
            title(currentField)
            if any(strcmp(currentField, {'spatialBinInfoRMSD'} ))
                colormap(ax, flipud(parula))
            end
        end
    else
        warning('Missing desired field in op for plot')
    end
end

%% Plot Place field scores
if analysisParam.calcPFscores && (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)

    analysisTitle = 'PF score (optimal=1, saturated<0)';
    cbLabel1 = 'Mean PF score';
    cbLabel2 = 'Best PF score';

    figure;
    imagesc(xParamvec, yParamvec, -1*mean([op.PFscore], 3)', 'AlphaData', ~isnan( mean([op.PFscore], 3)' )&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = cbLabel1;
    xlabel(xName,'Interpreter','none')
    ylabel(yName,'Interpreter','none')
    title(analysisTitle)
    lims = caxis;
    caxis([max(0, lims(1)), inf])

    figure;
    imagesc(xParamvec, yParamvec, -1*min([op.PFscore], [], 3)', 'AlphaData', ~isnan( min([op.PFscore], [], 3)' )&~excludedBadConnProb')
    set(gca,'YDir','normal')
    cb = colorbar(); cb.Label.String = cbLabel2;
    xlabel(xName,'Interpreter','none')
    ylabel(yName,'Interpreter','none')
    title(analysisTitle)
    lims = caxis;
    caxis([max(0, lims(1)), inf])
end


%% Plot place field stat. distributions for an example parameter point
% Plot all stats, without averaging over cells. Matches exptAnalysisPlaceFields.m

if analysisParam.calcPFscores && (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
    title_list = {'Peak rate', 'Specificity', 'Spatial Info.'};
    plotFields = {'Peak rate', 'Specificity', 'Spatial Info.'};
    myPlotSettings(width=6, height=2);
    figure;
    t = tiledlayout('flow', 'Padding', 'compact'); % 'TileSpacing','compact'
    title(t, ['Analysis over cells'])
    toPlot = [1:3];
    for ithStat = 1:numel(title_list)
        X_cells = squeeze(op.exampleCells(ithStat,:,:,:));
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
end


%% Plot place field stat. distributions for an example parameter point
% Plot all stats, means over PFmaps

if analysisParam.calcPFscores && (analysisParam.plotExtraFigs || ~analysisParam.MSPlots)
    myPlotSettings(width=7.0, height=3.5);

    title_list = {'Peak rate', 'Specificity', 'Spatial Info.', 'PF Score (optimal=1)', ...
        'Peak RMSD', 'Mean PF RMSD', 'Frac center 1/3rd', 'binSpatialInfo'};
    plotFields = {'meanPFPeaks', 'meanCellSpecificity', 'meanSpatialInfo', 'PFscore', ...
        'linearPeakRMSD', 'meanPFRMSD', 'fracPeakCenter3rd', 'spatialBinInfoMean'};

    figure;
    t = tiledlayout('flow', 'Padding', 'compact'); % 'TileSpacing','compact'
    title(t, ['Analysis over populations'])
    for ithStat = 1:numel(plotFields)
        %X = tempOP(ithStat,:);
        X = op.(plotFields{ithStat})(analysisParam.paramPoint(1),analysisParam.paramPoint(2),:);
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

    myPlotSettings
end


%% Functions:

function op = loopOPFields(simParam, analysisParam, opFields, PFresultsStruct, modelParam)
% Loops over parameters and networks to run calcOPFields()
for ithParam1 = 1:numel(simParam.variedParam(1).range)
    for ithParam2 = 1:numel(simParam.variedParam(2).range)

        % Run calculations
        temp_op = cell(numel(opFields), simParam.nNets); % needed for parfor usage
        if isempty(PFresultsStruct(ithParam1, ithParam2, 1).results)
            for ii = 1:numel(temp_op); temp_op{ii} = [nan]; end
        else
            if analysisParam.useParfor
                parfor ithNet = 1:simParam.nNets
                    temp_op(:,ithNet) = calcOPFields(ithParam1, ithParam2, ithNet, analysisParam, opFields, PFresultsStruct, modelParam);
                end
            else
                for ithNet = 1:simParam.nNets
                    temp_op(:,ithNet) = calcOPFields(ithParam1, ithParam2, ithNet, analysisParam, opFields, PFresultsStruct, modelParam);
                end
            end % End calculations
        end

        % Move temp_op data into the op struct
        % The i dimension of temp_op(i,:) is in the order of the fields of opFields
        for ithField = 1:numel(opFields)
            % Check special cases
            if strcmp(opFields{ithField}, 'exampleCells')
                if ithParam1==analysisParam.paramPoint(1) && ithParam2==analysisParam.paramPoint(2)
                    op.exampleCells =  [temp_op{strcmp('exampleCells', opFields), :}];
                end
            else % Typical case
                op.(opFields{ithField})(ithParam1, ithParam2, :) = [temp_op{ithField,:}];
            end
        end

    end% ithParam2 loop
end % ithParam1 loop
end

function op_net = calcOPFields(ithParam1, ithParam2, ithNet, analysisParam, opFields, PFresultsStruct, modelParam)
% Calculates all desired statistics for a given net and param point
%
%

% Set up for analysis:
op_net = cell(numel(opFields), 1);

% Get matrix of PFs
if iscell(PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields)
    PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields{analysisParam.ithEnv};
    E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
elseif ismatrix(PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields) && analysisParam.ithEnv==1
    PFmat = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields;
    E_indices = PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
else
    warning(['PF data error ', num2str(ithParam1), ' ', num2str(ithParam2), ' ', num2str(ithNet)])
    return
end

% Select only active, excitatory cells
PFinds = [ ismember(1:size(PFmat, 1), E_indices)' & (max(PFmat, [], 2)>=analysisParam.minPeakRate) ];
ePFmat = PFmat(PFinds,:);
iPFmat = PFmat(~ismember(1:size(PFmat, 1), E_indices)' ,:);
%{
figure; histogram(max(ePFmat, [], 2))
figure; histogram(ePFmat)
figure; histogram(max(iPFmat, [], 2))
figure; histogram(iPFmat)
%}
% Calculate PF spatial distribution statistics
nBins = size(ePFmat, 2); % number of spatial bins
nPCs = size(ePFmat, 1); % number of place cells

[~, pfPeakSpatialInd] = max(ePFmat, [], 2);
[peakBinCounts,EDGES] = histcounts(pfPeakSpatialInd, 0.5:1:nBins+0.5);

% fields
for ithField = 1:numel(opFields)
    currentField = opFields{ithField};

    switch currentField
        case 'meanPFPeaks'
            % 'MeanRate' of PF peak
            cellPeaks = max(ePFmat, [], 2);
            %fprintf('PF mean Peak Rate %0.4f %1.0f %1.0f %1.0f \n', meanPeakRate, ithParam1, ithParam2, ithNet)
            op_net{ithField} = mean(cellPeaks);

        case 'meanIRate'
            % 'MeanIRate'
            op_net{ithField} = mean(iPFmat, 'all');
            %{
            [v,i] = max(accum_ksstat); x = ( ePFmat(i,:)-mean(ePFmat(i,:), 2) )./std(ePFmat(i,:));
            figure; plot( ( ePFmat(i,:)-mean(ePFmat(i,:), 2) )./std(ePFmat(i,:)) ); title(['KS-stat ', num2str(v)])
            [F,X] = ecdf(x); figure; hold on; plot(X, cdf(makedist('Normal'), X)); ecdf(x); title(['KS-stat ', num2str(v)])
            [v,i] = min(accum_ksstat); x = ( ePFmat(i,:)-mean(ePFmat(i,:), 2) )./std(ePFmat(i,:));
            figure; plot( ( ePFmat(i,:)-mean(ePFmat(i,:), 2) )./std(ePFmat(i,:)) ); title(['KS-stat ', num2str(v)])
            [F,X] = ecdf(x); figure; hold on; plot(X, cdf(makedist('Normal'), X)); ecdf(x); title(['KS-stat ', num2str(v)])
            keyboard
            %}

        case 'meanCellSpecificity'
            % 'specificity'
            cellSpecificity = 1 - mean( ePFmat>[0.25*max(ePFmat, [], 2)], 2 );
            %fprintf('PF specificity %0.4f %1.0f %1.0f %1.0f \n', mean(cellSpecificity), ithParam1, ithParam2, ithNet)
            op_net{ithField} = mean(cellSpecificity);

        case 'meanSpatialInfo'
            % 'information'
            spatialInfo = nanmean( [ePFmat./mean(ePFmat, 2)] .* log(( ePFmat+eps )./mean(ePFmat, 2) ), 2 );
            %fprintf('PF info. %0.4f %1.0f %1.0f %1.0f \n', nanmean(spatialInfo), ithParam1, ithParam2, ithNet)
            op_net{ithField} = nanmean(spatialInfo);

        case 'linearPeakRMSD'
            % 'linear peaks rmsd'
            op_net{ithField} = sqrt(mean( (peakBinCounts-(nPCs/nBins)).^2 ) );


        case 'peakDistKLD'
            % 'linear peaks KL-divergence'
            p_data = peakBinCounts./sum(peakBinCounts);
            p_uniforma = ones(size(p_data)) * 1/numel(p_data);
            op_net{ithField} = -nanmean( p_data .* log2(p_uniforma./p_data) ) * numel(p_data);

        case 'meanPFRMSD'
            % 'Mean PF rmsd'
            popPF = mean(ePFmat, 1);
            op_net{ithField} = sqrt(mean( (popPF-(nPCs/nBins)).^2 ) );

        case 'fracPeakCenter3rd'
            % Frac peaks in center 1/3rd
            firstInd = round(nBins/3*1);
            secondInd = round(nBins/3*2);
            op_net{ithField} = sum(peakBinCounts(firstInd:secondInd))/nPCs;

        case 'spatialBinInfoMean'
            % 'rmsd of the spatial information in each spatial bin'
            % or 'bottom x%-tile of bin spatial info'
            % or 'variance of bin spatial info'
            % or 'mean bin spatial info', this correlates with cell-wise spatial info
            binSpatialInfo = nanmean( [ePFmat./mean(ePFmat, 1)] .* log(( ePFmat+eps )./mean(ePFmat, 1) ), 1 );
            op_net{ithField} = nanmean(binSpatialInfo);

        case 'spatialBinInfo33pct'
            binSpatialInfo = nanmean( [ePFmat./mean(ePFmat, 1)] .* log(( ePFmat+eps )./mean(ePFmat, 1) ), 1 );
            op_net{ithField} = prctile(binSpatialInfo, 33);

        case 'spatialBinInfoVar'
            binSpatialInfo = nanmean( [ePFmat./mean(ePFmat, 1)] .* log(( ePFmat+eps )./mean(ePFmat, 1) ), 1 );
            op_net{ithField} = var(binSpatialInfo);

        case 'spatialBinInfoRMSD'
            binSpatialInfo = nanmean( [ePFmat./mean(ePFmat, 1)] .* log(( ePFmat+eps )./mean(ePFmat, 1) ), 1 );
            op_net{ithField} = sqrt(mean( (binSpatialInfo-mean(binSpatialInfo)).^2 ) );

        case 'PFscore'
            % Calculate place field score
            PFscore=nan;
            if analysisParam.calcPFscores
                % Reformat matrix of PFs to struct required for calculate_linfieldsScore()
                network = struct;
                linfields = {};
                network.E_indices = E_indices;
                network.all_indices = 1:modelParam.n;
                day = 1; epoch = 1; tetrode = 1; tr = 1;
                linfields{day}{epoch}{tetrode}{1}{tr}(:,1) = modelParam.gridxvals*100; % convert to cm
                for ithCell = network.E_indices
                    linfields{day}{epoch}{tetrode}{ithCell}{tr}(:,5) = PFmat(ithCell,:);
                end
                PFscore = calculate_linfieldsScore(linfields, modelParam, network, trajectory=ithEnv);
                fprintf('PF score %0.4f %1.0f %1.0f %1.0f \n', PFscore, ithParam1, ithParam2, ithNet)
            end
            op_net{ithField}=PFscore;

        case 'exampleCells'
            if analysisParam.plotExampleParamPoint && ithParam1==analysisParam.paramPoint(1) && ithParam2==analysisParam.paramPoint(2)
                cellPeaks = max(ePFmat, [], 2);
                cellSpecificity = mean( ePFmat<=[0.25*max(ePFmat, [], 2)], 2 );
                spatialInfo = nanmean( [ePFmat./mean(ePFmat, 2)] .* log(( ePFmat+eps )./mean(ePFmat, 2) ), 2 );
                op_net{ithField} = cell(3,1);
                op_net{ithField}{1} = cellPeaks;
                op_net{ithField}{2} = cellSpecificity;
                op_net{ithField}{3} = spatialInfo;
            end

        otherwise
            error('Unkonwn fieldname for op in opFields')

    end % field switch
end % opFields loop
end % end function