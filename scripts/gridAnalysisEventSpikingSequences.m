%% Analyze and plot within-event spiking and decode statisics across parameter grids
% gridAnalysisEventSpiking.m
%

%% Choose which simulation results to analyze

decodeFileID = "mainGrid300s"; % Select which decoding file result ID to analyze
%decodeFileID = "mainPaperGrid"; % Select which decoding file result ID to analyze


%% Set up and load data

myPlotSettings

addpath(fullfile("..", "src"))
Config = utils.setConfig;
dataPath = Config.simDataPath;

% Load parameters and results
decodeName = grid.getDecodeFilename(decodeFileID);
if exist('pfGridspikes', 'var')
    disp('Assuming loaded variables')
else
    [modelParam, ...
        simParam, ...
        resultsStruct, ...
        PFresultsStruct, ...
        preplayGridSpikes, ...
        pfGridspikes ...
        ] = grid.loadResults(decodeName);
end


%% Analysis parameters
analysisParam.paramPoint = [2, 3]; % Select which example param point to plot, can be "all"
analysisParam.analyzeEntireGrid = false; % If false, only analyze analysisParam.paramPoint
analysisParam.ithEnv = 1; % which environment's decoding and place fields to plot/analyze
analysisParam.exampleEvent = [2, 3, 1, 1]; % Indices of example event to plot (param1, param2, net, event)

if analysisParam.analyzeEntireGrid
    analysisParam.paramSetInds = combvec(1:size(resultsStruct, 1), 1:size(resultsStruct, 2))';
else
    analysisParam.paramSetInds = analysisParam.paramPoint;
end

% Analysis parameters
analysisParam.minPfRate = 3; % Exclude cells with PF peaks below this rate
analysisParam.pfRankMethod = 'mean'; % peak or mean
analysisParam.invertDirectionality = false; % if true, invert time in events with negative decode slopes
analysisParam.onlySingleClustCells = false; % if true, exclude cells that belong to multiple clusters from the analysis

analysisParam.usePfCmUnits = true;

analysisParam.normEventLength = true;
analysisParam.smoothWindow = round(10 * 1e-3 / modelParam.dt); %round(1/10 * numel(event_timeVec));
analysisParam.activeClustMuliplyer = 2; % Clusters must be this times more active than any other cluster to be called 'active'
analysisParam.calcCorrForEachNet = false;

% Options to run/plot
analysisParam.plotExtraFigs = false;
analysisParam.runRasterChecks = false; % For each event, verify that the participating cells in the raster match the decode cells

analysisParam.figSettings = 'manuscript'; % 'manuscript' or ''

disp(analysisParam)

%% Parameter grid loop: analysis and plotting

switch analysisParam.figSettings
    case 'manuscript'
        myPlotSettings(width=2.0, height=1.75)
    otherwise
        myPlotSettings
end

tic
rng('default');
for ithParamSetToPlot = 1:size(analysisParam.paramSetInds, 1)

    % Set up varied-parameters for this loop
    param1Ind = analysisParam.paramSetInds(ithParamSetToPlot,1);
    param1Val = simParam.variedParam(1).range(param1Ind);
    disp(['Param1=', num2str(param1Ind), ...
        ': ', simParam.variedParam(1).name, '=', num2str(param1Val)])
    param2Ind = analysisParam.paramSetInds(ithParamSetToPlot,2);
    param2Val = simParam.variedParam(2).range(param2Ind);
    disp(['Param2=', num2str(param2Ind), ...
        ': ', simParam.variedParam(2).name, '=', num2str(param2Val)])
    linearParamInd = find( ...
        all([param1Val; param2Val] == simParam.parameterSets_vec, 1) ...
        );
    if isempty(linearParamInd)
        continue
    end

    % Update parameters for the ithParamSetToPlot:
    loopModelParam = modelParam;
    loopModelParam.(simParam.variedParam(1).name) = simParam.variedParam(1).range(param1Ind);
    loopModelParam.(simParam.variedParam(2).name) = simParam.variedParam(2).range(param2Ind);
    loopModelParam = set_depedent_parameters(loopModelParam);

    % Run analysis calculation
    output = grid.calculateEventSpikingSequences(...
        analysisParam, ...
        loopModelParam, ...
        resultsStruct(param1Ind, param2Ind, :), ...
        PFresultsStruct(param1Ind, param2Ind, :), ...
        preplayGridSpikes{linearParamInd} ...
        );

end

runTime = toc;
disp(['Runtime: ', datestr(datenum(0,0,0,0,0,runTime),'HH:MM:SS')])

