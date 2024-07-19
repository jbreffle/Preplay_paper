%% This script analyzes PF map correlations for PFs calculated from subsets of trials
% gridAnalysisPlaceFieldReliability.m
%
%

myPlotSettings

%% Set up and load data

% Choose results to analyze
decodeFileID = "mainPaperGrid"; % Select which decoding file result ID to analyze
% Param set [2,3] for (mcn=1.25, cluster=15)
% Param set [9,3] for (mcn=3.0, cluster=15)
analysisParam.ithSimParam = [2, 3]; % Select which example param point to plot

addpath(fullfile("..", "src"))
Config = utils.setConfig;
dataPath = Config.simDataPath;

% Load chosen data
analysisParam.simName = grid.getDecodeFilename(decodeFileID);
tmp = split(analysisParam.simName, 'decode');
simName = tmp{1};
load([dataPath filesep simName filesep simName 'parameters.mat'])
%load([dataPath filesep simName filesep analysisParam.simName '.mat'])
%load([dataPath filesep simName filesep simName 'spikes_preplay.mat'])
load([dataPath filesep simName filesep simName 'spikes_PF.mat'])

disp(['Analyzing ', analysisParam.simName])
param1Val = simParam.variedParam(1).range(analysisParam.ithSimParam(1));
disp(['Param1=', num2str(analysisParam.ithSimParam(1)), ': ', simParam.variedParam(1).name, '=', num2str(param1Val)])
param2Val = simParam.variedParam(2).range(analysisParam.ithSimParam(2));
disp(['Param2=', num2str(analysisParam.ithSimParam(2)), ': ', simParam.variedParam(2).name, '=', num2str(param2Val)])

linearParamInd = find(all([param1Val; param2Val] == simParam.parameterSets_vec, 1));

PFSpikeTimes_paramPoint = PFSpikeTimes_grid{linearParamInd};
rngStruct_paramPoint = rngStruct_grid{linearParamInd};
clear PFSpikeTimes_grid rngStruct_grid


%% Analysis parameters

analysisParam.maxNSplits = 6; % Maximum number of indepednent trial splits (other than just odd/even)
analysisParam.useSplitNForAllTrajCorrs = true;

analysisParam.figSettings = 'standard'; % standard, manuscript

analysisParam.saveResults = false;
analysisParam.plotResults = true;


%% Perform PF correlation analyses, for each network

% Initialize the main output structures of the analysis
allTrajCorrs = nan(simParam.nNets, modelParam.nEnvironments, modelParam.nEnvironments); % All pair-wise trajectory PF map comparisons
allHalfCorrs = nan(simParam.nNets, modelParam.nEnvironments); % odd/even PF map correlations
allDirCorrs = nan(simParam.nNets, modelParam.nEnvironments); % left/right PF map correlations

nSplitTrials = floor(modelParam.nTrials_PF/2); % number of trials used for each split
nPossibleSplits = nchoosek(2*nSplitTrials, nSplitTrials);
nActualSplits = min(nPossibleSplits, analysisParam.maxNSplits);
allSplitHalfCorrs = nan(simParam.nNets, modelParam.nEnvironments*nActualSplits); % odd/even PF map correlations

tic
for ithNet = 1:simParam.nNets

    %% Set up network struct
    network.E_indices = rngStruct_paramPoint.net(ithNet).E_indices;

    %% Convert spike times to binary spike matrix
    netSpikeTimes = reshape(PFSpikeTimes_paramPoint(ithNet,:,:,:), [modelParam.nEnvironments, modelParam.nTrials_PF, modelParam.n]);
    %opS = false(pfsim.nEnvironments, pfsim.nTrials, parameters.n, pfsim.t_steps);
    opS = false(modelParam.n, numel(modelParam.t_PF), modelParam.nEnvironments, modelParam.nTrials_PF);
    for ithEnv = 1:modelParam.nEnvironments
        for ithTrial = 1:modelParam.nTrials_PF
            for ithCell = 1:modelParam.n %E_indices
                % cellSpikeInds = round(netSpikeTimes{ithCell}./parameters.dt);
                cellSpikeInds = round(netSpikeTimes{ithEnv, ithTrial, ithCell}./modelParam.dt);
                opS(ithCell, cellSpikeInds, ithEnv, ithTrial) = 1;
            end
        end
    end

    %% Place Field correlation analysis

    % Calculate odd/even trial linfields
    oddTrialInds = 1:2:(2*nSplitTrials);
    evenTrialInds = 2:2:(2*nSplitTrials);

    linfieldsOdd = calculate_linfields(opS(:,:,:,oddTrialInds), modelParam);
    [PFmatOdd, PFpeaksSequenceOdd] = calculate_PFmat(linfieldsOdd, modelParam, network, plotPFsFlag=false);

    linfieldsEven = calculate_linfields(opS(:,:,:,evenTrialInds), modelParam);
    [PFmatEven, PFpeaksSequenceEven] = calculate_PFmat(linfieldsEven, modelParam, network, plotPFsFlag=false);

    % left/right correlations, for both odd and even trial PFs
    tmp = [];
    for ithEnv = 2:2:modelParam.nEnvironments
        % Note: PFmat includes I-cell rows, and need to reverse the PFs for opposite directions
        pixelCorrsOdd = diag(corr(PFmatOdd{ithEnv-1}(network.E_indices,:), fliplr(PFmatOdd{ithEnv}(network.E_indices,:)) ));
        mapCorrOdd = nanmean(pixelCorrsOdd);

        pixelCorrsEven = diag(corr(PFmatEven{ithEnv-1}(network.E_indices,:),  fliplr(PFmatEven{ithEnv}(network.E_indices,:)) ));
        mapCorrEven = nanmean(pixelCorrsEven);

        tmp = [tmp, mapCorrOdd, mapCorrEven];
    end
    allDirCorrs(ithNet,:) = tmp;


    % odd/even correlations, for each traj
    tmp = [];
    for ithEnv = 1:1:modelParam.nEnvironments
        pixelCorrs = diag(corr(PFmatOdd{ithEnv}(network.E_indices,:), PFmatEven{ithEnv}(network.E_indices,:) ));
        mapCorr = nanmean(pixelCorrs);
        tmp = [tmp, mapCorr];
    end
    allHalfCorrs(ithNet,:) = tmp;


    % odd/even correlations, for each traj and all possible trial divisions
    tmp = [];
    allSplitCombs = utils.splitVector(1:(2*nSplitTrials));
    for ithSplit = 1:nActualSplits
        % Calculate trial-split linfields
        split1Inds = allSplitCombs{ithSplit}{1};
        split2Inds = allSplitCombs{ithSplit}{2};

        linfields1 = calculate_linfields(opS(:,:,:,split1Inds), modelParam);
        [PFmat1, PFpeaksSequence1] = calculate_PFmat(linfields1, modelParam, network, plotPFsFlag=false);

        linfields2 = calculate_linfields(opS(:,:,:,split2Inds), modelParam);
        [PFmat2, PFpeaksSequence2] = calculate_PFmat(linfields2, modelParam, network, plotPFsFlag=false);

        for ithEnv = 1:modelParam.nEnvironments
            pixelCorrs = diag(corr(PFmat1{ithEnv}(network.E_indices,:), PFmat2{ithEnv}(network.E_indices,:) ));
            mapCorr = nanmean(pixelCorrs);
            tmp = [tmp, mapCorr];
        end

    end
    allSplitHalfCorrs(ithNet,:) = tmp;



    % Calculate all within and across env corrs
    if analysisParam.useSplitNForAllTrajCorrs
        linfields = calculate_linfields(opS(:,:,:,1:nSplitTrials), modelParam);
    else
        linfields = calculate_linfields(opS, modelParam);
    end
    [PFmat, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network, plotPFsFlag=false);
    for ithEnv = 1:(modelParam.nEnvironments-1)

        ithEnvPFmap = PFmat{ithEnv}(network.E_indices,:);
        if mod(ithEnv, 2)==1
            ithEnvPFmap = fliplr(ithEnvPFmap);
        end

        for jthEnv = (ithEnv+1):(modelParam.nEnvironments)

            jthEnvPFmap = PFmat{jthEnv}(network.E_indices,:);
            if mod(jthEnv, 2)==1
                jthEnvPFmap = fliplr(jthEnvPFmap);
            end

            pixelCorrs = diag(corr(ithEnvPFmap, jthEnvPFmap));
            mapCorr = nanmean(pixelCorrs);
            allTrajCorrs(ithNet,ithEnv,jthEnv) = mapCorr;
        end
    end

end

runTime = toc;
disp(['Runtime: ', num2str(runTime/60), ' min'])


%% Save results

if analysisParam.saveResults
    %{
    fileName = [simName, 'variedDecode', datestr(now,'yyyy-mm-ddTHH-MM')];
    save_path = fullfile(dataPath, simName, fileName);
    save(save_path, 'PFresultsStruct_variedDecode', 'resultsStruct_variedDecode', 'analysisParam', 'modelParam')
    %}
    disp('DID NOT SAVE RESULTS, need to write the code')
else
    disp('DID NOT SAVE RESULTS')
end


%% Plot results

if analysisParam.plotResults

    switch analysisParam.figSettings
        case 'standard'
            myPlotSettings(width=3.5, height=2.5)
        case 'manuscript'
            myPlotSettings(width=2.5, height=1.5)
    end

    allCorr = [allHalfCorrs(:); allDirCorrs(:)];

    % analysisParam.nBins = 2*ceil(sqrt(numel([allHalfCorrs(:); allDirCorrs(:)])));
    analysisParam.nHistBins = 15;

    %figure; histogram(allHalfCorrs, nBins); ylabel('All odd/even comp. (count)'); xlabel('PF map correlation')
    %figure; histogram(allSplitHalfCorrs, nBins); ylabel('All odd/even comp. (count)'); xlabel('PF map correlation')
    %figure; histogram(allDirCorrs, nBins); ylabel('All left/right comp. (count)'); xlabel('PF map correlation')
    %figure; histogram(allCorr, nBins); ylabel('All comp. (count)'); xlabel('PF map correlation')

    % Compare odd/even and left/right
    [~,EDGES] = histcounts(allCorr, analysisParam.nHistBins);
    figure; hold on
    histogram(allHalfCorrs, EDGES)
    histogram(allDirCorrs, EDGES)
    legend({'Odd/Even','Left/Right'}, 'location', 'north', 'Box', 'off')
    ylabel('All pair-wise (count)'); xlabel('Remapping correlation (r)')

    [~, pval, ~] = kstest2(allHalfCorrs(:), allDirCorrs(:));
    figure; hold on
    ecdf(allHalfCorrs(:))
    ecdf(allDirCorrs(:))
    title(['kstest pval=', num2str(pval)])
    legend({'Odd/Even','Left/Right'}, 'location', 'best', 'Box', 'off')
    ylabel('CDF'); xlabel('Remapping correlation (r)')


    % Compare all possible splits and left/right
    [~,EDGES2] = histcounts([allSplitHalfCorrs(:); allDirCorrs(:)], analysisParam.nHistBins);
    figure; hold on
    histogram(allSplitHalfCorrs, EDGES2, normalization='probability')
    histogram(allDirCorrs, EDGES2, normalization='probability')
    legend({'All splits','Left/Right'}, 'location', 'best', 'Box', 'off')
    ylabel('All pair-wise (prob.)'); xlabel('Remapping correlation (r)')

    [~, pval, ~] = kstest2(allSplitHalfCorrs(:), allDirCorrs(:));
    figure; hold on
    ecdf(allSplitHalfCorrs(:))
    ecdf(allDirCorrs(:))
    title(['kstest pval=', num2str(pval)])
    legend({'All splits','Left/Right'}, 'location', 'best', 'Box', 'off')
    ylabel('CDF'); xlabel('Remapping correlation (r)')


    % Compare all pair-wise traj using complete trials for all PFs
    [~,EDGES] = histcounts([allTrajCorrs(:); allSplitHalfCorrs(:)], analysisParam.nHistBins);
    figure; hold on
    histogram(allTrajCorrs(~isnan(allTrajCorrs)), EDGES, normalization='probability')
    histogram(allSplitHalfCorrs, EDGES, normalization='probability')
    if analysisParam.useSplitNForAllTrajCorrs
        legend({'All traj (matched nTrials)', 'Odd/Even trial PFs'}, location='north', Box='off')
    else
        legend({'All traj (all nTrials)', 'Odd/Even trial PFs'}, location='north', Box='off')
    end
    ylabel('All pair-wise (prob.)'); xlabel('Remapping correlation (r)')

end

