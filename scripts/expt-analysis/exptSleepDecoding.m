%% Code to run sleep-session decoding of Shin et al., 2019 data
% sleep_decode_batch_new.m
%
% This was adapted from exptSleepDecoding.m to run decoding of
% experimental data from Shin et al., 2019 using the decode_events function
% created for simulated preplay events.


%% Set up

% Get data directory
userDir = char(java.lang.System.getProperty('user.home'));
if ispc; boxDir = 'Box';
elseif ismac; boxDir = ['Library' filesep 'CloudStorage' filesep 'Box-Box'];
else; disp('error');
end

addpath(fullfile("..", "..", "src"))
Config = utils.setConfig;
% TODO: switch to using Config struct

animalDataDir = 'Data/ReplayNet-Data/js_SingleDayExpt';

myPlotSettings


%% Parameters needed for the decode_events() function, to match analysis of simulated data

% Decoding parameters, taken from simulation analysis
modelParam.cellcountthresh = 5;     % Minimum number of participating cells for a candidate event
modelParam.shuffleIterations = 100;	% Number of time-bin shuffles for each event
modelParam.useLogSumDecode = true;	% Use sum of log probabilities for decoding events
modelParam.tBinSz = 10;            % ms, time bin size for decoding
modelParam.minEventDur = 50;       % ms, exclude events shorter than this
modelParam.wellcutoff = 0;         % cm, remove reward-well regions (15cm around start and end); or 0 cm without exclusion
modelParam.minPeakRate = 3;        % Hz, minimum peak rate to include cell as Place Cell
modelParam.normByTraj_decode = true; % Should be true (false is equivalent to old method, which was wrong)
modelParam.downSampleCells = false; % true/false flag to use downsampling
modelParam.downSampledFrac = 0.1;  % fraction of cells to use in decode, if downSampleCells==true
modelParam.dispFlag_decode = 1; % display event analysis progression if ==1

% new options
modelParam.useMESrippletime = true; % Use the *rippletimes_MES.mat files for decoding


%% Set up for analysis

animalprefix_list = {'ER1','KL8','JS14','JS15','JS17','JS21'};
% animalprefix_list = {'JS17'};
day_list = 1;           % Vector of days to process
%eps = 2:2:16;           % Vector of epochs to analyze
eps = [1,2,3];
%eps = 1;           % Vector of epochs to analyze

% Directory to load data from to run replay_decoding_CA1_singleday()
loaddir = fullfile(userDir, boxDir, animalDataDir);

% Directory to save results of replay_decoding_CA1_singleday()
savedir = fullfile(userDir, boxDir, animalDataDir, 'Preplay_decoding', filesep, datestr(now,'yyyy-mm-ddTHH-MM'), filesep);

% Parameters/options
savedata = true;           % save data = 1; not save = 0
figopt = 1;             % 1 = plot decoding result for each event; 0 = no figure generated


%% Loop over animals, days, and epochs

% Create the directory for saving decoded event data
if ~exist(savedir) && savedata
    mkdir(savedir)
end

tic
mainScriptSeed = 0; % rng(0)==rng('default')
rng(mainScriptSeed)
for ithAnimal = 1:length(animalprefix_list)
    animalprefix = animalprefix_list{ithAnimal};
    for day = day_list
        for epoch = eps

            % replay_decoding_CA1_singleday(animalprefix,day,ep,cellcountthresh,loaddir,savedir,savedata,figopt)

            decode_events(modelParam, animalprefix, ...
                savedir, savedata, day=day, epoch=epoch, figopt=figopt);

            disp(['Running time: ', num2str(toc/60/60), ' hours'])
        end
    end
end
runTime = toc;

batch_decode_git_hash_string = system('git rev-parse HEAD');
batch_decode_filename = mfilename('fullpath');

if savedata
    save([savedir, 'mainParams.mat'])
end