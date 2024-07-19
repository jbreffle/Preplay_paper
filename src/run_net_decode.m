function replaytrajectory = run_net_decode(modelParam, network, trialResults, linfields, preplaySpikeMat, varargin)
% run_net_decode.m
%
% Runs preplay decoding on data from one simulated network
%
% Inputs:
%   - modelParam: the model parameter structure
%   - network: the network structure
%   - trialResults: output of detect_PBE()
%   - linfields: Jadhav lab formated linear PF struct
%   - preplaySpikeMat: Binary spike matrix (can be either all cells or just
%                           E-cells)
%
% Outputs:
%   - replaytrajectory: Bayesian decoding results cell array.
%       - replaytrajectory{p.day}{p.ep}.<fieldname>
%
% Optional:
%   -
%
% Side effects:
%   - Creates then deletes temporary files in ./results/tmp/
%   - If decodeFuncSeed is given, rng is modified
%   - If saveDecode is give, decoding results are saved
%   - If fitopt is give, figures are plotted
%

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'modelParam',     @isstruct)
addRequired(inputObj, 'network',        @isstruct)
addRequired(inputObj, 'trialResults',	@isstruct)
addRequired(inputObj, 'linfields',      @iscell)
addRequired(inputObj, 'preplaySpikeMat', @islogical)
addParameter(inputObj, 'day',   1, @isnumeric);
addParameter(inputObj, 'ep',    2, @isnumeric);
addParameter(inputObj, 'eprun', 2, @isnumeric);
addParameter(inputObj, 'tet',   1, @isnumeric);
addParameter(inputObj, 'figopt',     0, @isnumeric);  % 2=new plot for each decode, 1=plot in single figure
addParameter(inputObj, 'saveDecode', false, @islogical);
addParameter(inputObj, 'decodeFuncSeed', -1, @isnumeric);
parse(inputObj, modelParam, network, trialResults, linfields, preplaySpikeMat, varargin{:});
p = inputObj.Results;


%% Set up data structures needed for Jadhav lab code, from randnet_PF.m sim

tetinfo{p.day}{p.eprun}{p.tet}.numcells = modelParam.n_E;

if numel(linfields{p.day})<p.eprun
    linfields{p.day}{p.eprun} = linfields{p.day}{1};
end
if numel(linfields{p.day}{p.eprun}{:})==modelParam.n
    linfields{p.day}{p.eprun}{p.tet}(network.I_indices) = [];
end

ripple{p.day}{p.eprun}.starttime = trialResults.events(:,1).*modelParam.dt;
ripple{p.day}{p.eprun}.endtime = trialResults.events(:,2).*modelParam.dt;

if size(preplaySpikeMat, 1)==numel(network.E_indices)
    E_Inds = 1:numel(network.E_indices);
else
    E_Inds = network.E_indices;
end

for cellid = 1:numel(network.E_indices)
    spikes{p.day}{p.eprun}{p.tet}{cellid}.data(:,1) = find(preplaySpikeMat(E_Inds(cellid),:)).*modelParam.dt; % spike times
end

%% Save SJlab formated files to tmp dir

[filepath,name] = fileparts( tempname([pwd, '/results/tmp/']) );
tmpSimPrefix = name; save_path = [filepath, filesep, tmpSimPrefix, '_direct'];
mkdir(save_path)

save([filepath, filesep, tmpSimPrefix, '_direct/', tmpSimPrefix, 'tetinfo.mat'], 'tetinfo')
save([filepath, filesep, tmpSimPrefix, '_direct/', tmpSimPrefix, 'linfields01.mat'], 'linfields')
save([filepath, filesep, tmpSimPrefix, '_direct/', tmpSimPrefix, 'rippletime01.mat'], 'ripple')
save([filepath, filesep, tmpSimPrefix, '_direct/', tmpSimPrefix, 'spikes01.mat'], 'spikes')


%% Run preplay decoding analysis

warning off
replaytrajectory = decode_events(modelParam, ...
    tmpSimPrefix, ...
    [filepath, filesep], ...
    p.saveDecode, ...
    day=p.day, ...
    epoch=p.ep, ...
    figopt=p.figopt, ...
    decodeFuncSeed=p.decodeFuncSeed, ...
    network=network ...
    );
warning on


%% Delete tmp dir files:
rmdir(save_path, 's')


end