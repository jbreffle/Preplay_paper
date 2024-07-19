%% Compares the spatial information vs the temporal information from simulated cells
% 1) Simulate the same networks with different durations of place field
% trajectory traversals
% 2) Combine the activity from the different trials by either time or
% space, then calculate the resulting information


%% Set up

addpath(fullfile("../src"))
Config = utils.setConfig;
baseDir = Config.simDataPath;

% Select which trials speed simulations to compare
comparisonDurs = [2, 4];

% Define the filenames: Identical networks, 2s, 1s, and 4s long place field trajectories
pf4sFile = 'grid_2024-06-10T12-26';
pf2sFile = 'grid_2024-06-04T07-47';

% Analysis/plotting options
ithParam1 = 1;
ithParam2 = 1;
ithTraj = 1;
plotSettings = 'manuscript';
plotExtra = false;
markExampleCell = false;

% Exclude cells whose average firing rate was below this during the PF simulations
minMeanRate = 0.5;

calcFieldInfo = @(X) nanmean( [X./nanmean(X, 2)] .* log(( X+eps )./nanmean(X, 2) ), 2 );

switch plotSettings
    case 'manuscript'
        myPlotSettings(width=3, height=2)
    otherwise
        myPlotSettings
end


%% Set up and load

assert(all(ismember(comparisonDurs, [2, 4])))

simADur = comparisonDurs(1);
simBDur = comparisonDurs(2);

switch comparisonDurs(1)
    case 2
        pfAFile = pf2sFile;
    case 4
        pfAFile = pf4sFile;
end

switch comparisonDurs(2)
    case 2
        pfBFile = pf2sFile;
    case 4
        pfBFile = pf4sFile;
end

resultsA = load(fullfile(baseDir, pfAFile, [pfAFile, 'results.mat']));
resultsB = load(fullfile(baseDir, pfBFile, [pfBFile, 'results.mat']));

paramsA = load(fullfile(baseDir, pfAFile, [pfAFile, 'parameters.mat']));
paramsB = load(fullfile(baseDir, pfBFile, [pfBFile, 'parameters.mat']));

assert(isequal(paramsA.simParam.t_max_PF, simADur))
assert(isequal(paramsB.simParam.t_max_PF, simBDur))

assert(isequal(paramsA.modelParam.t_max_PF, simADur))
assert(isequal(paramsB.modelParam.t_max_PF, simBDur))


%% Run analysis, looping over networks

nNets = size(resultsB.PFresultsStruct, 3);

placeInfo = [];
timeInfo = [];

allPlaceFieldsMat = [];
allTimeFieldsMat = [];

allMats.pfA = [];
allMats.pfB = [];
allMats.tfA = [];
allMats.tfB = [];
allMats.pfComb = [];
allMats.tfComb = [];

for ithNet = 1:nNets

    eIndsB = resultsB.PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
    eIndsA = resultsA.PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.E_indices;
    assert(isequal(eIndsA, eIndsA)) % Networks should be identical

    pfMatB = resultsB.PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields{ithTraj};
    pfMatA = resultsA.PFresultsStruct(ithParam1, ithParam2, ithNet).results{1}.linfields{ithTraj};

    nCells = size(pfMatB, 1);
    nBins = size(pfMatB, 2);

    %% Pad to equivalent lengths
    % For short simulations take mean of adjacent columns and combine, then pad with nans

    tfMatB = pfMatB;
    tfMatA = nan(size(pfMatA));
    for ithBin = 1:2:nBins
        tfMatA(:,ceil(ithBin/2)) = (pfMatA(:,ithBin)+pfMatA(:,ithBin+1))/2;
    end

    %% Calculate average firing rate: mean of linfield (weighted by run duration if combining)

    meanRateB = nanmean(tfMatB, 2);
    meanRateA = nanmean(tfMatA, 2);
    meanRateAll = (meanRateB*simBDur + meanRateA*simADur) / (simBDur+simADur);

    %% Combine into all-trial place fields (PFs) and all-trial time fields (TFs)

    validIndsB = ismember(eIndsB, find(meanRateAll>minMeanRate));
    validIndsA = ismember(eIndsA, find(meanRateAll>minMeanRate));
    assert(isequal(validIndsA, validIndsB))

    meanPFMat = (pfMatB(validIndsB,:)  + pfMatA(validIndsA,:) ) / 2;
    meanTimeFieldMat = (tfMatB(validIndsB,:) + tfMatA(validIndsA,:) ) / 2;

    placeInfo = [placeInfo; calcFieldInfo(meanPFMat)];
    timeInfo = [timeInfo; calcFieldInfo(meanTimeFieldMat)];

    allPlaceFieldsMat = [allPlaceFieldsMat; meanPFMat];
    allTimeFieldsMat = [allTimeFieldsMat; meanTimeFieldMat];

    allMats.pfA = [allMats.pfA; pfMatA(validIndsA,:)];
    allMats.pfB = [allMats.pfB; pfMatB(validIndsB,:)];
    allMats.tfA = [allMats.tfA; tfMatA(validIndsA,:)];
    allMats.tfB = [allMats.tfB; tfMatB(validIndsB,:)];
    allMats.pfComb = [allMats.pfComb; meanPFMat];
    allMats.tfComb = [allMats.tfComb; meanTimeFieldMat];

    if 0
        cellInd = 5;
        figure; hold on
        plot(tfMatA(cellInd,:))
        plot(tfMatB(cellInd,:))
        figure; hold on
        plot(pfMatA(cellInd,:))
        plot(pfMatB(cellInd,:))
    end

end


%% Plot place and time fields

if plotExtra
    % figure; imagesc(allPlaceFieldsMat)

    plotNormRates = true;
    posBins = 1:2:100;
    nCells = size(allPlaceFieldsMat, 1);

    % For each trajectory, plot the PFs and the distribution of peaks
    row_all_zeros = find(all( allPlaceFieldsMat==0, 2)) ;
    row_n_all_zeros = find(~all( allPlaceFieldsMat==0, 2)) ;
    [~,peakRateLocation] = max(squeeze(allPlaceFieldsMat(row_n_all_zeros,:)), [], 2);
    [~,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    curr_PFpeaksSequence = [row_n_all_zeros(sortedCellIndsbyPeakRateLocation); row_all_zeros];
    peakRates = max(allPlaceFieldsMat(curr_PFpeaksSequence,:), [], 2);
    if plotNormRates
        rateDenom1 = peakRates;
        caxmax = 1;
    else
        rateDenom1 = 1;
        caxmax = max(allPlaceFieldsMat, [], 'all');
    end
    figure; tiledlayout(1,2); sgtitle(['Env ID ', num2str(ithTraj)], fontweight='normal', fontsize=12)
    nexttile
    imagesc(posBins, 1:nCells, allPlaceFieldsMat(curr_PFpeaksSequence,:)./rateDenom1 );
    colorbar;
    clim([0, caxmax])
    xlabel('Position (cm)'); ylabel('Cell (sorted)');
    nexttile
    histogram(peakRates)
    xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');


    % figure; imagesc(allTimeFieldsMat)

    % For each trajectory, plot the PFs and the distribution of peaks
    row_all_zeros = find(all( allTimeFieldsMat==0, 2)) ;
    row_n_all_zeros = find(~all( allTimeFieldsMat==0, 2)) ;
    [~,peakRateLocation] = max(squeeze(allTimeFieldsMat(row_n_all_zeros,:)), [], 2);
    [~,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    curr_PFpeaksSequence = [row_n_all_zeros(sortedCellIndsbyPeakRateLocation); row_all_zeros];
    peakRates = max(allTimeFieldsMat(curr_PFpeaksSequence,:), [], 2);
    if plotNormRates
        rateDenom1 = peakRates;
        caxmax = 1;
    else
        rateDenom1 = 1;
        caxmax = max(allTimeFieldsMat, [], 'all');
    end
    figure; tiledlayout(1,2); sgtitle(['Env ID ', num2str(ithTraj)], fontweight='normal', fontsize=12)
    nexttile
    imagesc(posBins, 1:nCells, allTimeFieldsMat(curr_PFpeaksSequence,:)./rateDenom1 );
    colorbar;
    clim([0, caxmax])
    xlabel('Position (cm)'); ylabel('Cell (sorted)');
    nexttile
    histogram(peakRates)
    xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');
end


%% Plot spatial and temporal information

myPlotSettings(width=1.75, height=1.5)

medianPlaceInfo = nanmedian(placeInfo)
medianTimeInfo = nanmedian(timeInfo)
size(placeInfo)

figure(634);
hold on;
scatter(placeInfo, timeInfo, 10)
ref1 = refline(1);
ref1.Color = 'k';
ref1.LineStyle = '--';
% ref2 = refline(0.5); ref2.Color = 'k'; ref2.LineStyle = '--';
scatter(medianPlaceInfo, medianTimeInfo, 'filled', 'ro')
xlabel('Place information')
ylabel({' ', 'Time information'})
title(' ')
xlim([0, max([timeInfo; placeInfo])*1.05])
ylim([0, max([timeInfo; placeInfo])*1.05])


figure
hold on
ecdf(placeInfo)
ecdf(timeInfo)
xlabel('Information (bits/spike)')
ylabel({'Cumulative', 'proportion'})
cdfLegend = legend({'Place', 'Time'}, 'location', 'best');
cdfLegend.ItemTokenSize = [30*0.5, 18*1.0];
[H,P,KSSTAT] = kstest2(placeInfo, timeInfo);
title(['p=', num2str(P, 2)]) %, 'Interpreter', 'latex')


%% Plot example cell's time and place fields

myPlotSettings(width=3.25, height=1.5)

if plotExtra
    infoDiff = placeInfo-timeInfo;
    [val, idx] = max(infoDiff);
    [sortVals, sortInds] = sort(infoDiff, 'descend');

    figure;
    hold on;
    plot(allPlaceFieldsMat(idx,:))
    plot(allTimeFieldsMat(idx,:))

    %{
    figure; hold on
    plot(allMats.pfA(idx,:))
    plot(allMats.pfB(idx,:))
    
    figure; hold on
    plot(allMats.tfA(idx,:))
    plot(allMats.tfB(idx,:))
    
    figure; hold on
    plot(allMats.tfA(sortInds(1000),:))
    plot(allMats.tfB(sortInds(1000),:))
    %}
end

nSteps = 50;
xTime = simBDur/nSteps:simBDur/nSteps:simBDur;
xPlace = 2:2:100;

% Cell with median place info:
if isequal(minMeanRate, 0.5)
    idx = 1361; % median cell when minMeanRate=2
else
    [y, idx] = min(abs(placeInfo-nanmedian(placeInfo)));
end
%idx = 1100

placeInfo(idx)
timeInfo(idx)
figure;
%
nexttile
hold on
plot(xPlace, allMats.pfA(idx,:), 'k-')
plot(xPlace, allMats.pfB(idx,:), 'k:')
ylabel('Firing rate (Hz)')
xlabel('Location (cm)')
%
nexttile
hold on
plot(xTime, allMats.tfA(idx,:), 'k-')
plot(xTime, allMats.tfB(idx,:), 'k:')
xlabel('Time (s)')
%
%nexttile; hold on
%plot(allMats.pfComb(idx,:))
%plot(allMats.tfComb(idx,:))
%
exampleLegend = legend({'2 s', '4 s'}, 'location', 'best');
exampleLegend.ItemTokenSize = [30*0.5, 18*1.0];
sgtitle(['Place=', num2str(placeInfo(idx)), ' time=', num2str(timeInfo(idx))], ...
    'fontsize', 10)


%% Mark example cell on the scatter plot

if markExampleCell
    % Values for cell with median placeInfo when minMeanRate=2.0
    placeInfoTarget = 0.48746;
    timeInfoTarget = 0.52726;
    [val1, idx2] = min(abs(placeInfo-placeInfoTarget));
    [val2, idx3] = min(abs(timeInfo-timeInfoTarget));
    assert(isequal(idx2, idx3));

    figure(634)
    hold on
    scatter(placeInfo(idx2), timeInfo(idx3), 'k', 'filled')
end
