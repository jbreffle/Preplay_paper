function replaytrajectory = decode_events(modelParam, animalprefix, savedir, savedata, varargin)
% Adapted from code for Shin et al., 2019, Jadhav Lab
% Previous function name was preplay_decoding_CA1_singleday.m
% Original function name from Wenbo was replay_decoding_CA1_singleday.m
%
% INPUTS:
%    animalprefix = animal prefix.
%    day = experimental day.
%    ep = epoch.
%    cellcountthresh = mininum number of cells active for considering as a
%                      cadidate event, usually = 5
%    savedir = directory for saving results
%    savedata = save results, 1 = save, 0 = not save
%    figopt = plot figures, 1 = plot, 0 = no plot

% This function requires the files:
% cellinfo, linfields01, rippletime01, spikes01, and tetinfo
%
% Producing rippletime01 requires the script replay_SWRtime_singleday.m,
% which requires the files:
% pos01, ripples01, tetinfo
%
% Assumptions:
% - all simulation data is on tetrode 1
% - all simulation data does not exclude any cells
%
% Loads: spikes, tetinfo, linfields, and rippletime


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'animalprefix',	@ischar)
addRequired(inputObj, 'savedir',	@ischar)
addRequired(inputObj, 'savedata',	@islogical)
addParameter(inputObj, 'day', 1,	@isnumeric)
addParameter(inputObj, 'epoch',	2, @isnumeric)
addParameter(inputObj, 'fix_maxJump', true,	@islogical)  % Old default was implicitly false
addParameter(inputObj, 'figopt', 0,	@isnumeric)  % 2=new plot for each decode, 1=plot in single figure
addParameter(inputObj, 'decodeFuncSeed',  -1,  @isnumeric)
addParameter(inputObj, 'network',  struct,  @isstruct)
parse(inputObj, modelParam, animalprefix, savedir, savedata, varargin{:});
p = inputObj.Results;

day = p.day;
epoch = p.epoch;
figopt = p.figopt;

% Only set the rng if a seed was passed in
if ~isequal(p.decodeFuncSeed, -1)
    rng(p.decodeFuncSeed, 'twister');
end

% Get analysis parameters from modelParam struct
tBinSz = modelParam.tBinSz; % ms, time bin size for decoding
minEventDur = modelParam.minEventDur; % ms, exclude events shorter than this
wellcutoff = modelParam.wellcutoff; % cm, remove reward-well regions (15cm around start and end); or 0cm without exclusion
minPeakRate = modelParam.minPeakRate; % Hz, minimum peak rate to include cell as Place Cell
downSampleCells = modelParam.downSampleCells; % true/false flag to use downsampling
downSampledFrac = modelParam.downSampledFrac; % fraction of cells to use in decode, if downSampleCells==1
dispFlag_decode = modelParam.dispFlag_decode; % display event analysis progression if ==1
shuffleIterations = modelParam.shuffleIterations; % Number of time-bin shuffles for each event
cellcountthresh = modelParam.cellcountthresh; % Minimum number of participating cells for a candidate event

msToSConversion = 1 / 1000;
sToMsConversion = 1000;

if isfield(modelParam, 'normByTraj_decode') && ~modelParam.normByTraj_decode
    warning("decode_events: modelParam.normByTraj_decode=false is deprecated ")
end
if isfield(modelParam, 'useLogSumDecode') && ~modelParam.useLogSumDecode
    warning("decode_events: modelParam.useLogSumDecode=false is deprecated ")
end


%% Set epochs
if ~mod(epoch,2)
    eprun = epoch; % using current epoch to construct place-field templates for RUN epoch
elseif epoch == 1
    eprun = 2; % using the first run epoch for pre-task sleep epoch
else
    eprun = epoch - 1; % using the previous run epoch for sleep epoch (post-task sleep)
end


%% Setup load/save paths

singleDayAnimals = {'ER1','KL8','JS14','JS15','JS17','JS21'};
socialWAnimals = {'XFB3'};
allAnimals = [singleDayAnimals, socialWAnimals];

isSimulatedData = ~ismember(animalprefix, allAnimals);

if isSimulatedData
    loadDir = fullfile(savedir, [animalprefix, '_direct/']);
    saveDir = loadDir;
else
    savedir_parts = strsplit(savedir, filesep);
    loadDir = fullfile(savedir_parts{1:end-3},  [animalprefix, '_direct/']);
    saveDir = savedir;
end


%% Get hpidx
% hpidx(ithCell, :) = [tetInd, cellInd]

[~, hpidx] = decode.getCellIdx( ...
    loadDir, animalprefix, day, epoch, isSimulatedData ...
    );

% Make sure only E-cells were included in the linfields struct, if sim data
if isSimulatedData
    assert(size(hpidx,1)==modelParam.n_E)
end

% If downSampleCells, use only random subset of cells for decoding
if downSampleCells
    nEpochCells = size(hpidx, 1);
    nSubsampleCells = round((1-downSampledFrac) * nEpochCells);
    randSubsetInds = sort(randsample(nEpochCells, nSubsampleCells));
    hpidx = hpidx(randSubsetInds, :);
end


%% Final set up

% Create the spatial ratemaps
[rm, pm, tm, cellidxm] = decode.createSpatialMatrices( ...
    animalprefix, hpidx, loadDir, day, eprun, modelParam ...
    );

% Load spikes
spikes = loaddatastruct(loadDir, animalprefix, 'spikes', day);

% Get riptimes
riptimes = decode.getRipTimes( ...
    animalprefix, loadDir, day, epoch, modelParam ...
    );

% Constants
hpnum = length(rm(1,:)); % update cell number (cell with peak rate < minPeakRate excluded)
nTraj = numel(unique(tm)); % N trajectories
nPBin = size(rm,1); % N positional bin
trajinfo = mean(tm,2); % trajactory number
expectedSpikes = rm .* (tBinSz * msToSConversion); % [nPos x nCells] Expected number of spikes per time bin
expectedSpikes  = reshape(expectedSpikes, [nPBin,1, hpnum]); %[nPos x 1 x nCell]
exponExpectedSpikes = exp(-expectedSpikes); %Exponent of equation.


%% Cell shuffle, if selected

% Note: if decodeCellShuffle==None, shuffledCellidxm = cellidxm
if isfield(modelParam, 'decodeCellShuffle')
    shuffledCellidxm = decode.shuffleCellidxm(cellidxm, modelParam, p.network);
    unShuffledSpikes = spikes;
    for ithCell = 1:hpnum
        spikes{day}{epoch}{cellidxm(ithCell,1)}{cellidxm(ithCell,2)} = ...
            unShuffledSpikes{day}{epoch}{shuffledCellidxm(ithCell,1)}{shuffledCellidxm(ithCell,2)};
    end
    clear unShuffledSpikes
end


%% Get cell data

[celldata, eventindex] = decode.getCellEventData( ...
    day, epoch, modelParam, riptimes, cellidxm, spikes ...
    );


%% Decoding loop

for event = 1:length(eventindex)

    % Show current event number
    if dispFlag_decode
        disp(['Event ', num2str(event), ' of ', num2str(length(eventindex))])
    end

    % Spikes per cell per time bin for current event
    spkPerBin = decode.getBinnedEventSpikes( ...
        celldata, riptimes, eventindex, event, modelParam ...
        );
    nSpkPerTBin = squeeze(sum(spkPerBin,3)); %[nTBin x 1] n spikes in tBin
    nonzerobins = find(nSpkPerTBin > 0);

    if figopt==1
        myPlotSettings(width=nTraj*2, height=2)
        figure(9342) % Re-plot new event in same figure
    elseif figopt==2
        myPlotSettings(width=nTraj*2, height=2)
        figure; % Plot in new figure
    end

    % trajectory loop (all trajectories stored in linfields{day}{epoch}{tetrode}{cell}{traj})
    for traj = 1:nTraj

        trajidx = find(trajinfo == traj);
        distvector = pm(trajidx,1)';

        % Calculate posterior probability matrix
        pMat = decode.calculatePosteriorMatrix( ...
            expectedSpikes(trajidx,:,:), exponExpectedSpikes(trajidx,:,:), spkPerBin ...
            );

        szPM2 = size(pMat,2);

        % Monte Carlo simulation for linear regression
        rvalues = [];
        slopes = [];
        entropy = [];
        interc = [];
        fitpVals = [];
        totalsamples = 10000;
        for rloop = 1:500
            tBinPicks = distsample(totalsamples,nSpkPerTBin);
            regressdata = [];
            entropy_loop = zeros(1, length(nonzerobins));
            for i = 1:length(nonzerobins)
                if (nSpkPerTBin(nonzerobins(i)) > 0)
                    tmpnumsamples = sum(tBinPicks == nonzerobins(i));
                    if ~isempty(find(pMat(:,nonzerobins(i)) ~= 0))
                        distpicks = distvector(distsample(tmpnumsamples,pMat(:,nonzerobins(i))))';
                        histBinNorm = hist(distpicks,0:5:200)./length(distpicks);
                        entropy_loop(i) = -nansum(histBinNorm.*log(histBinNorm));
                        distpicks(:,2) = nonzerobins(i);
                        regressdata = [regressdata; distpicks];
                    end
                end
            end
            regressdata(:,3) = 1;

            [b,bint,r,rint,stats] = regress(regressdata(:,1),[regressdata(:,3),regressdata(:,2)]);

            rvalues = [rvalues; stats(1)];
            fitpVals = [fitpVals; stats(3)];
            slopes = [slopes; b(2)];
            interc = [interc; b(1)];
            entropy = [entropy; mean(entropy_loop)];
        end
        % WIP: move above into a function, return a struct?
        %{
        figure; plot(regressdata)
        figure; scatter(regressdata(:,2), regressdata(:,1))
        figure; histogram(rvalues)
        keyboard
        %}

        % Calculate max jump distance
        [~, PeakDecodeInd] = max(pMat, [], 1);
        if p.fix_maxJump
            % Remove bins with no spikes for calculating maxJump,
            % unless fix_maxJump=false was given as an optional input
            PeakDecodeInd(sum(pMat, 1)==0) = [];
        end
        if ~isempty(diff(PeakDecodeInd))
            maxJump(event,traj) = max(abs(diff(PeakDecodeInd))./size(pMat, 1));
        else
            maxJump(event,traj) = nan;
        end

        % Store results
        Res(event,traj) =  mean(rvalues);
        Fit_pVal(event, traj) = mean(fitpVals);
        YInterc(event,traj) = mean(interc);
        Spd(event,traj) = mean(slopes);
        Entropy(event,traj) = mean(entropy);
        pMat_cell{event}{traj}.pMat = pMat;
        pMat_cell{event}{traj}.timevec = 1:szPM2;
        pMat_cell{event}{traj}.posvec = distvector;
        pMat_cell{event}{traj}.timebinsz = tBinSz;

        % Calculate weighted linear regression
        [X,y] = meshgrid([1:size(pMat,2)],[1:size(pMat,1)]);
        w = pMat;
        mdl = fitlm(X(:),y(:),'Weights',w(:));
        weightedSlope(event,traj) = mdl.Coefficients.Estimate(2);
        weightedR2(event,traj) = mdl.Rsquared.Ordinary;

        if figopt
            subplot(1,nTraj,traj)
            imagesc(1:szPM2,distvector,pMat);
            colormap(jet)
            hold on
            plot([1, szPM2], [0, Spd(event,traj) * (szPM2-1)] + (YInterc(event,traj) +1), 'w','linewidth',2);
            hold off
            %sgtitle([animalprefix, 'Day: ', num2str(day), ', Ep: ', num2str(epoch), ', Event:', num2str(event)]);
            sgtitle(['Event: ', num2str(event) , ' of ' , num2str(length(eventindex))]);
            caxis([0 0.1])
        end

        %-------Shuffling to get the pvalue for each traj------%
        permbins = nonzerobins;
        srvalues = [];
        smaxJumps = [];
        sslopes = [];
        sweightedSlope = [];
        sweightedR2 = [];
        for iteration = 1:shuffleIterations % 1500
            permbins = permbins(randperm(length(permbins)));% temporal shuffle
            % calculate shuffled pMat
            tmpspkPerBin = zeros(size(spkPerBin));
            tmpspkPerBin(:,permbins,:) = spkPerBin(:,nonzerobins,:);
            tmpfactSpkPerBin = factorial(tmpspkPerBin); %Factorial to divide by
            tmpnSpkPerTBin = squeeze(sum(tmpspkPerBin,3)); %[nTBin x 1] number of spikes in tBin

            tmpexpecSpk = expectedSpikes(trajidx,:,:);
            tmpexpon = exponExpectedSpikes(trajidx,:,:);
            wrking_shuff = bsxfun(@power, tmpexpecSpk, tmpspkPerBin); %[nPos x nTbin x nCell]
            wrking_shuff = bsxfun(@rdivide, wrking_shuff, tmpfactSpkPerBin); %[nPos x nTbin x nCell]
            wrking_shuff = bsxfun(@times,wrking_shuff, tmpexpon); %[nPos x nTbin x nCell]
            logTmppMat = sum( log(wrking_shuff) , 3); % log transformed, non-normalized prob [nPos x Tbin]
            tmppMat = exp(logTmppMat-max(logTmppMat)); %Peak-normalized prob [nPos x Tbin]
            tmppMat(:,tmpnSpkPerTBin==0)  =0; % so the posterior matrix can be smoothed.
            tmppMat(isnan(tmppMat)) = 0;
            for i = 1:szPM2 % normalized across positions to 1 for each time bin
                if (sum(tmppMat(:,i))>0)
                    tmppMat(:,i) = tmppMat(:,i)./sum(tmppMat(:,i));
                end
                if (sum(tmppMat(:,i))==0)
                    %disp('Shuffle event with a zero row')
                end
            end
            clear wrking_shuff tmpfactSpkPerBin tmpexpon

            % Monte Carlo simulation for linear regression
            tBinPicks = distsample(totalsamples,nSpkPerTBin);
            regressdata = [];
            for i = 1:length(permbins)
                if (nSpkPerTBin(nonzerobins(i)) > 0)
                    tmpnumsamples = sum(tBinPicks == nonzerobins(i));
                    if ~isempty(find(pMat(:,nonzerobins(i)) ~= 0))
                        distpicks = distvector(distsample(tmpnumsamples,tmppMat(:,permbins(i))))';
                        distpicks(:,2) = permbins(i);
                        regressdata = [regressdata; distpicks];
                    end
                end
            end
            regressdata(:,3) = 1;
            [b,bint,r,rint,stats] = regress(regressdata(:,1),[regressdata(:,3),regressdata(:,2)]);
            srvalues = [srvalues; stats(1)];
            sslopes = [sslopes; b(2)];

            [~, PeakDecodeInd] = max(tmppMat, [], 1);
            if p.fix_maxJump
                % Remove bins with no spikes for calculating maxJump,
                % unless fix_maxJump=false was given as an optional input
                PeakDecodeInd(sum(tmppMat, 1)==0) = [];
            end
            currentShuffMaxJump = max(abs(diff(PeakDecodeInd))./size(tmppMat, 1));
            smaxJumps = [smaxJumps, currentShuffMaxJump];

            [X,y] = meshgrid([1:size(tmppMat,2)],[1:size(tmppMat,1)]);
            w = tmppMat;
            mdl = fitlm(X(:),y(:),'Weights',w(:));
            sweightedSlope = [sweightedSlope,  mdl.Coefficients.Estimate(2)];
            sweightedR2 = [sweightedR2, mdl.Rsquared.Ordinary];

        end

        % calculate p-value
        pvalue(event,traj) = sum(Res(event,traj) < srvalues)/length(srvalues);

        % Save shuffle stats
        shuffle_Spd{event}{traj} = sslopes;
        shuffle_rvalues{event}{traj} = srvalues;
        shuffle_maxJump{event}{traj} = smaxJumps;
        shuffle_weightedSlope{event}{traj} = sweightedSlope;
        shuffle_weightedR2{event}{traj} = sweightedR2;
    end

    [minP,tidx] = min(pvalue(event,:));% find the minimun pvalue
    if minP < 0.05 % significant
        decode_traj(event) = tidx; % trajectory with the minimun pvalue represented
    else
        decode_traj(event) = 0;% no significant traj
    end

    % Save cell info during the event
    cellsiAll = celldata(find(celldata(:,2)==eventindex(event)),3);
    [cellsi, ~] = unique(cellsiAll, 'first');
    activecell{event} = cellsi;
    activecellidx{event} = cellidxm(cellsi,:);

    if figopt==2
        for subplotind = 1:nTraj
            subplot(1,nTraj,subplotind);
            xlabel('Time (10 ms bin)');
            title(['p-val=', num2str(pvalue(event, subplotind))]);
            if subplotind == 1
                ylabel('Position (cm)');
            end
        end
        drawnow
    end

end


%% Store output struct data

% structure result
if ~exist('pMat_cell')
    % Original output
    replaytraj.pMat         = [];
    replaytraj.maxJump      = [];
    replaytraj.rsquare      = [];
    replaytraj.FitpVal      = [];
    replaytraj.slopes       = [];
    replaytraj.YInterc      = [];
    replaytraj.Entropy      = [];
    replaytraj.eventidx     = [];
    replaytraj.besttraj     = [];
    replaytraj.pvalue       = [];
    replaytraj.activecell   = [];
    replaytraj.activecellidx = [];
    replaytraj.sigeventprc  = [];
    replaytraj.sigeventnum  = [];
    replaytraj.candeventnum = [];
    % Shuffle data
    replaytraj.shuffle_rsquare = [];
    replaytraj.shuffle_maxJump = [];
    replaytraj.shuffle_slopes  = [];
    % New, slopes and weighted decodes
    replaytraj.weightedSlope    = [];
    replaytraj.weightedR2       = [];
    replaytraj.shuffle_weightedSlope = [];
    replaytraj.shuffle_weightedR2   = [];
    % New,
    replaytraj.eventStart = [];
    replaytraj.eventEnd = [];
else
    % Original output
    replaytraj.pMat         = pMat_cell;
    replaytraj.maxJump      = maxJump; % maximum jump in peak decoded pos. between time bins, in frac of track
    replaytraj.rsquare      = Res;
    replaytraj.FitpVal      = Fit_pVal;
    replaytraj.slopes       = Spd;
    replaytraj.YInterc      = YInterc;
    replaytraj.Entropy      = Entropy;
    replaytraj.eventidx     = eventindex;
    replaytraj.besttraj     = decode_traj;
    replaytraj.pvalue       = pvalue;
    replaytraj.activecell   = activecell;
    replaytraj.activecellidx = activecellidx;
    replaytraj.sigeventprc  = length(find(decode_traj~=0))./length(decode_traj);
    replaytraj.sigeventnum  = length(find(decode_traj~=0));
    replaytraj.candeventnum = length(decode_traj);
    % Shuffle data
    replaytraj.shuffle_rsquare = shuffle_rvalues;
    replaytraj.shuffle_maxJump = shuffle_maxJump;
    replaytraj.shuffle_slopes   = shuffle_Spd;
    % New, weighted decode correlation
    replaytraj.weightedSlope    = weightedSlope;
    replaytraj.weightedR2       = weightedR2;
    replaytraj.shuffle_weightedSlope	= shuffle_weightedSlope;
    replaytraj.shuffle_weightedR2       = shuffle_weightedR2;
    % New,
    replaytraj.eventStart = riptimes(:,1);
    replaytraj.eventEnd = riptimes(:,2);
end
replaytraj.wellcutoff = wellcutoff;
replaytraj.tBinSz = tBinSz;
replaytraj.cellcountthresh = cellcountthresh;
replaytraj.decodeFuncSeed = p.decodeFuncSeed;

if isfield(modelParam, 'decodeCellShuffle')
    replaytraj.decodeCellShuffle = modelParam.decodeCellShuffle;
    replaytraj.shuffledCellidxm = shuffledCellidxm;
    replaytraj.cellidxm = cellidxm;
end


replaytrajectory{day}{epoch} = replaytraj;


%% Save data

if savedata
    save(fullfile(sprintf('%s%sreplaydecode_CA1_%02d_%02d.mat', saveDir,animalprefix,day,epoch)), 'replaytrajectory');
end


%% Return to original plotting settings

if figopt
    myPlotSettings()
end

