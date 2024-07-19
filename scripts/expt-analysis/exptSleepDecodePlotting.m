%% Sleep decode analysis of Shin et al., 2019
% exptSleepDecodePlotting.m
%
% Statistical analysis of the decoded events from the first sleep session
% of Shin et al., 2019


%% Set up

% Get data directory
userDir = char(java.lang.System.getProperty('user.home'));
if ispc; boxDir = 'Box';
elseif ismac; boxDir = ['Library' filesep 'CloudStorage' filesep 'Box-Box'];
else; disp('error');
end

addpath(fullfile("..", "..", "src"))
Config = utils.setConfig;

%animalDataDir = 'Data/js-single-day-data/Preplay_decoding/';
animalDataDir = Config.exptPreplayDecodePath;

myPlotSettings(width=2.75, height=2)


%% analysis set up

animal_names = Config.animalPrefixes;
%animal_names = {'ER1','KL8','JS14','JS15','JS17','JS21'};

decodeSet{1} = '2023-09-19T16-16'; % fixed rng, epochs 1:3, decode_events(), MES rippletimes, after fixing maxJump

% Options
decodeSetToPlot = 1;
analysisParam.trajToInclude = [1:4];
% analysisParam.trajToInclude = 1;

animalprefix_list = strcat(decodeSet{decodeSetToPlot}, filesep, animal_names);
analysisParam.animalprefix_list = animalprefix_list;
analysisParam.useWeightedCorr = true;
analysisParam.ithDay = 1;
analysisParam.ithEpoch = 1;
analysisParam.plotIndividualAnimals = false;
analysisParam.timeBinWidth = 10/1000; % 10 ms

analysisParam.bestEventMethod = 'r'; % pval, r, jd
analysisParam.nBestEvents = 24; % Size adjusts for 12, 24, 48, or 96
analysisParam.minActiveBins = 5; % Don't plot example events with fewer than this active bin count

analysisParam.onlyBestTraj = false;
if analysisParam.onlyBestTraj; analysisParam.trajToInclude= [1:4]; end

analysisParam.plotPvalMat = true;
analysisParam.figSettings = 'manuscript'; % 'manuscript', 'manuscriptAlt' (for supplement)

r2_actual_all = [];
r2_shuffle_all = [];
jd_actual_all = [];
jd_shuffle_all = [];

bestEvents.pval = inf(analysisParam.nBestEvents, 1);
bestEvents.r = -inf(analysisParam.nBestEvents, 1);
bestEvents.jd = inf(analysisParam.nBestEvents, 1);
bestEvents.pMat = cell(analysisParam.nBestEvents, 1);

event_durations = [];
nEvents = zeros(size(animalprefix_list));

% mainParamStruct = load(fullfile(... 'mainParams.mat'));
disp(['Plotting ', decodeSet{decodeSetToPlot}, ', epoch ', num2str(analysisParam.ithEpoch), ', traj ', num2str(analysisParam.trajToInclude)])

%% Loop over animals

tic
for ithAnimal = 1:numel(animalprefix_list)

    % Load data
    fileLocation = fullfile(animalDataDir, ...
        [animalprefix_list{ithAnimal}, 'replaydecode_CA1_', num2str(analysisParam.ithDay, '%02.f'), '_', num2str(analysisParam.ithEpoch, '%02.f'), '.mat']);
    load(fileLocation)

    % Get actual and shuffle r2 for all trajectories
    if analysisParam.useWeightedCorr
        r2_actual_animal = replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.weightedR2;
        r2_shuffle_animal = vertcat( replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.shuffle_weightedR2{:});
    else
        r2_actual_animal = replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.rsquare;
        r2_shuffle_animal = vertcat( replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.shuffle_rsquare{:});
    end

    % Get actual and shufle jd for all trajectories
    jd_actual_animal = replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.maxJump;
    jd_shuffle_animal = vertcat( replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.shuffle_maxJump{:} );

    % Calculate each events pval of its r2 relative to its own shuffles
    pval_actual_animal = nan(size(r2_actual_animal));
    for ithEvent=1:size(r2_actual_animal, 1 )
        pval_actual_animal(ithEvent,:) = mean(r2_actual_animal(ithEvent,:)< vertcat(r2_shuffle_animal{ithEvent,:})', 1);
    end

    % Accumulate all of the best events
    switch analysisParam.bestEventMethod
        case 'jd'
            [sorted_pval, sorted_pval_ind] = sort(jd_actual_animal(:));
            for ithPval = 1:numel(sorted_pval)
                [tmpIndEvent, tmpIndTraj] = ind2sub(size(jd_actual_animal),sorted_pval_ind(ithPval));
                if any(sorted_pval(ithPval)<bestEvents.pval)
                    if sum(replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat{tmpIndEvent}{1}.pMat, 'all')>=analysisParam.minActiveBins
                        [~, tmpInd] = max(bestEvents.pval);
                        bestEvents.pMat{tmpInd} = replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat{tmpIndEvent}{tmpIndTraj}.pMat;
                        bestEvents.r(tmpInd) = sqrt(r2_actual_animal(tmpIndEvent, tmpIndTraj));
                        bestEvents.jd(tmpInd) = jd_actual_animal(tmpIndEvent, tmpIndTraj);
                        bestEvents.pval(tmpInd) = pval_actual_animal(tmpIndEvent, tmpIndTraj);
                    end
                else
                    break
                end
            end


        case 'r'
            [sorted_r, sorted_r_ind] = sort(sqrt(r2_actual_animal(:)), 'descend');
            for ithR = 1:numel(sorted_r)
                [tmpIndEvent, tmpIndTraj] = ind2sub(size(r2_actual_animal),sorted_r_ind(ithR));
                if any(sorted_r(ithR)>bestEvents.r)
                    if sum(replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat{tmpIndEvent}{1}.pMat, 'all')>=analysisParam.minActiveBins
                        [~, tmpInd] = min(bestEvents.r);
                        bestEvents.pMat{tmpInd} = replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat{tmpIndEvent}{tmpIndTraj}.pMat;
                        bestEvents.r(tmpInd) = sqrt(r2_actual_animal(tmpIndEvent, tmpIndTraj));
                        bestEvents.jd(tmpInd) = jd_actual_animal(tmpIndEvent, tmpIndTraj);
                        bestEvents.pval(tmpInd) = pval_actual_animal(tmpIndEvent, tmpIndTraj);
                    end
                else
                    break
                end
            end

            %[bestStat_animal, maxr2_linInd] = max(r2_actual_animal, [], 'all', 'linear');
            %[bestEventInd, bestTrajInd] = ind2sub(size(r2_actual_animal),maxr2_linInd);
        case 'pval'
            [sorted_pval, sorted_pval_ind] = sort(pval_actual_animal(:));
            for ithPval = 1:numel(sorted_pval)
                [tmpIndEvent, tmpIndTraj] = ind2sub(size(pval_actual_animal),sorted_pval_ind(ithPval));
                if any(sorted_pval(ithPval)<bestEvents.pval)
                    if sum(replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat{tmpIndEvent}{1}.pMat, 'all')>=analysisParam.minActiveBins
                        [~, tmpInd] = max(bestEvents.pval);
                        bestEvents.pMat{tmpInd} = replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat{tmpIndEvent}{tmpIndTraj}.pMat;
                        bestEvents.r(tmpInd) = sqrt(r2_actual_animal(tmpIndEvent, tmpIndTraj));
                        bestEvents.jd(tmpInd) = jd_actual_animal(tmpIndEvent, tmpIndTraj);
                        bestEvents.pval(tmpInd) = pval_actual_animal(tmpIndEvent, tmpIndTraj);
                    end
                else
                    break
                end
            end

            %[bestStat_animal, maxr2_linInd] = min(pval_actual_animal, [], 'all', 'linear');
            %[bestEventInd, bestTrajInd] = ind2sub(size(pval_actual_animal),maxr2_linInd);
        otherwise
            error('Unknown option for analysisParam.bestEventMethod')
    end
    % bestPmat_animal = replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat{bestEventInd}{bestTrajInd}.pMat;
    % figure; imagesc(bestPmat_animal); title(['Best decode for ', animal_names{ithAnimal}])

    % Select r2 and jd for only desired trajectories
    if analysisParam.onlyBestTraj
        % For each decoded event, select only the best trajectory stats
        % Actual event r2s
        [actual_min_r2s ,actual_min_r2_ind] = max(r2_actual_animal, [], 2);
        r2_actual_animal = actual_min_r2s;
        % Actual event jds
        tmpInd=sub2ind(size(jd_actual_animal),(1:size(jd_actual_animal,1))',actual_min_r2_ind);
        jd_actual_animal = jd_actual_animal(tmpInd);
        % [animal_actual_jumps, actual_min_r2_ind, animal_actual_jumps(tmpInd)] % verify with this line
        % Shuffle event r2s
        %tmp_r2 = cell(size(animal_actual_r2s));
        tmp_r2 = [];
        tmp_jd = [];
        for ithEvent = 1:size(r2_actual_animal, 1)
            [tmp_min, tmp_ind] = max(vertcat(r2_shuffle_animal{ithEvent,:})', [], 2);
            tmp_r2 = [tmp_r2; {tmp_min}];
            %tmp_r2{ithEvent} = {tmp_min};
            tmp_event_jd = vertcat(jd_shuffle_animal{ithEvent,:})';
            tmpInd=sub2ind(size(tmp_event_jd),(1:size(tmp_event_jd,1))',tmp_ind);
            tmp_jd = [tmp_jd; {tmp_event_jd(tmpInd)}];
            % [tmp_event_jd, tmp_ind, tmp_event_jd(tmpInd)] % verify with this line
        end
        r2_shuffle_animal = tmp_r2;
        jd_shuffle_animal = tmp_jd;
        % Shuffle event jds

    else
        % Take all trajectory states
        % Actual event r2s
        r2_actual_animal = r2_actual_animal(:,analysisParam.trajToInclude);
        r2_actual_animal = r2_actual_animal(:);
        % Actual event jds
        jd_actual_animal = jd_actual_animal(:,analysisParam.trajToInclude);
        jd_actual_animal = jd_actual_animal(:);
        % Shuffle event r2s
        r2_shuffle_animal = r2_shuffle_animal(:,analysisParam.trajToInclude);
        r2_shuffle_animal = r2_shuffle_animal(:);
        r2_shuffle_animal = cellfun(@transpose,r2_shuffle_animal,'UniformOutput',false);
        % Shuffle event jds
        jd_shuffle_animal = jd_shuffle_animal(:,analysisParam.trajToInclude);
        jd_shuffle_animal = jd_shuffle_animal(:);
        jd_shuffle_animal = cellfun(@transpose,jd_shuffle_animal,'UniformOutput',false);
    end

    % Store this animal's data
    nEvents(ithAnimal) = numel(r2_actual_animal);
    trajI = 1; % All traj event lengths are identical
    animal_event_durations = analysisParam.timeBinWidth.* cellfun(@(x) x{trajI}.timevec(end), replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.pMat);
    event_durations = [event_durations, animal_event_durations];

    r2_actual_all = [r2_actual_all; r2_actual_animal];
    r2_shuffle_all = [r2_shuffle_all; r2_shuffle_animal];

    jd_actual_all = [jd_actual_all; jd_actual_animal];
    jd_shuffle_all = [jd_shuffle_all; jd_shuffle_animal];

    % Calculate kstest for this animal
    [~,P,ksstat] = kstest2(sqrt(r2_actual_animal), vertcat(sqrt(vertcat(r2_shuffle_animal{:}))));

    % Plot this animal's data
    if analysisParam.plotIndividualAnimals
        figure; hold on;
        ecdf(sqrt(r2_actual_animal))
        ecdf(sqrt(vertcat(r2_shuffle_animal{:})))
        xlabel('|correlation|')
        ylabel('CDF')
        if analysisParam.onlyBestTraj
            title({[animal_names{ithAnimal}, ', epoch ', num2str(analysisParam.ithEpoch), ', best traj.'], ...
                [num2str(nEvents(ithAnimal)), ' events, pval: ', num2str(P)]}, ...
                fontWeight='normal')
        else
            title({[animal_names{ithAnimal}, ', epoch ', num2str(analysisParam.ithEpoch), ', traj. ', num2str(analysisParam.trajToInclude)], ...
                [num2str(nEvents(ithAnimal)), ' events, pval: ', num2str(P)]}, ...
                fontWeight='normal')
        end
    end

    disp(['Animal % sig. = ', num2str(replaytrajectory{analysisParam.ithDay}{analysisParam.ithEpoch}.sigeventprc)])
    disp(['Animal median duration = ', num2str( median(animal_event_durations))])
    disp(['Animal KS-test p-val = ', num2str(P)])
end

figure; hold on;
ecdf(sqrt(r2_actual_all))
ecdf(sqrt(vertcat(r2_shuffle_all{:})))
xlabel('|correlation|')
ylabel('CDF')
[~,P,ksstat] = kstest2(sqrt(r2_actual_all), sqrt(vertcat(r2_shuffle_all{:})));
if analysisParam.onlyBestTraj
    title({['All animals', ', epoch ', num2str(analysisParam.ithEpoch), ', best traj.'], ...
        [num2str(sum(nEvents(:))), ' events, pval: ', num2str(P)]}, ...
        fontWeight='normal')
else
    title({['All animals', ', epoch ', num2str(analysisParam.ithEpoch), ', traj. ', num2str(analysisParam.trajToInclude)], ...
        [num2str(sum(nEvents(:))), ' events, pval: ', num2str(P)]}, ...
        fontWeight='normal')
end
figure; histogram(event_durations);
title('All events')
xlabel('Ripple dur. (s)')
ylabel('Count (events)')

%{
figure; scatter(actual_jumps, sqrt(actual_r2s))
xlabel('jd'); ylabel('|r|'); fitlm(actual_jumps, sqrt(actual_r2s))
%}

%% Plot best events

tBinSz = 0.02; % 2 cm
spatialBin = 1;
cAxisRange = [0, 0.1];

if analysisParam.nBestEvents>0
    if analysisParam.nBestEvents==12
        myPlotSettings(width=3,height=3)
    elseif analysisParam.nBestEvents==24
        myPlotSettings(width=6,height=3)
    elseif analysisParam.nBestEvents==48
        myPlotSettings(width=6,height=6)
    elseif analysisParam.nBestEvents==96
        myPlotSettings(width=12,height=6)
    end

    switch analysisParam.bestEventMethod
        case 'r'
            [~, eventSortInd] = sort(bestEvents.r, 'descend');
        case 'pval'
            [~, eventSortInd] = sort(bestEvents.pval);
        case 'jd'
            [~, eventSortInd] = sort(bestEvents.jd);
    end

    figure; tfig = tiledlayout('flow', 'Padding', 'none', 'TileSpacing', 'compact');
    title(tfig, ['Best ', num2str(numel(bestEvents.pMat)), ' of ', num2str(sum(nEvents)), ' events'])
    for ithEventToPlot = eventSortInd'
        nexttile

        [yt, xt] = size(bestEvents.pMat{ithEventToPlot});
        imagesc([1:xt]*(tBinSz), [1:yt]*(spatialBin*100), bestEvents.pMat{ithEventToPlot})
        colormap hot
        title(['r=', num2str(bestEvents.r(ithEventToPlot), 2), '; jd=', num2str(bestEvents.jd(ithEventToPlot), 2)], 'FontWeight', 'normal')
        if isequal(analysisParam.bestEventMethod, 'pval')
            %title(['r=', num2str(bestEvents.r(ithEventToPlot), 2), '; jd=', num2str(bestEvents.jd(ithEventToPlot), 2), '; p=', num2str(bestEvents.pval(ithEventToPlot), 2) ], 'FontWeight', 'normal')
        end
        % caxis(([0, 0.75*max(bestEvents.pMat{ithEventToPlot}, [], 'all')]))
        caxis(cAxisRange)

        set(gca,'ytick',[])
        %yt = yticks; yticklabels(yt*(pfsim.spatialBin*100)); ylabel('Position (cm)')
        %xt = xticks; xticklabels(xt*(tBinSz)); xlabel('Time (ms)')
        %xticks('auto')
        %if find(ithEventToPlot==eventSortInd)<=(2*analysisParam.nEventsToPlot/3)
        %    set(gca,'xtick',[])
        %end
        set(gca,'xtick',[])
    end
    title(tfig, '')

    if ~exist('bestDecodeColorBarFigHandle', 'var')
        myPlotSettings(width=2, height=1.25)
        bestDecodeColorBarFigHandle = figure;
        title('colorbar', 'fontSize', 19); colorbar; colormap(hot); caxis(cAxisRange)
    end
end


%% Manuscript figures

op.allEventR2s = r2_actual_all;
op.allShuffleR2s = r2_shuffle_all; %  cellfun(@sqrt, op.allShuffleR2s, 'UniformOutput', false);
op.allEventMaxJumps = jd_actual_all;
op.allShuffleMaxJumps = jd_shuffle_all;

% Plot ECDF of r values
switch analysisParam.figSettings
    case 'standard'
        myPlotSettings
    case 'manuscript'
        myPlotSettings(width=2.5,height=1.5)
    case 'manuscriptAlt'
        myPlotSettings(width=1.6, height=1.25)
    case 'SfNPoster'
        myPlotSettings(width=6, height=3, lw=3, fzs=24, alw=3) % SfN-poster format
end
if ~isempty(op.allEventR2s)
    allEventR2s_ecdf = op.allEventR2s;
    r2vals_shuff_ecdf= (vertcat(op.allShuffleR2s{:}));

    figure; hold on; cdfplot(sqrt(allEventR2s_ecdf)); cdfplot(real(sqrt(r2vals_shuff_ecdf)));
    legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
    xlabel('|Correlation|'); ylabel({'Cumulative','proportion'});
    [~,P,KSSTAT] = kstest2(sqrt(abs(allEventR2s_ecdf)), sqrt(abs(r2vals_shuff_ecdf)) );
    disp(['Decode pval: ', num2str(P)])
    disp(['Decode ks-stat: ', num2str(KSSTAT)])
    disp(['Median shift: ', num2str(median(sqrt(allEventR2s_ecdf))-median(sqrt(r2vals_shuff_ecdf)))])
    title(['Experiments ', ...
        ' pval=', num2str(P), ...
        ' nEvents=', num2str(numel(allEventR2s_ecdf))])
    title ''
    legend({'Preplay', 'Shuffle'}, 'Location', 'Best'); legend boxoff
    grid off

    h = get(gca,'children');
    set(h(1), 'LineWidth', 1, 'color', [1, 0.3, 0.3], 'LineStyle', '-.')
    set(h(2), 'LineWidth', 1, 'color', [0, 0, 0.6])
    %set(h(2), 'LineWidth', 1, 'color', 'k', 'LineStyle', '--'); set(h(3), 'LineWidth', 1, 'color', [0, 0, 0.6])
else
    disp('No events to plot r-val CDF')
end

% Plot p-value matrix
if analysisParam.plotPvalMat
    switch analysisParam.figSettings
        case 'standard'
            myPlotSettings
        case 'manuscript'
            myPlotSettings(width=2.5, height=1.5)
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
        eventJD, shuffleJD, plotExtra=false, nanThreshIs0=true);

end

totalRunTime = toc
