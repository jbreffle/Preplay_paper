function linfields = calculate_linfields(opS, modelParam, varargin)
% calculate_linfields.m
%
% Calculate linear place fields, based on Jadhav lab methods.
% Adapted from ReplayNet code, with assumption of simulating a single
% trajectory for a single epoch.
%
% Inputs:
%   - opS: logical spike matrix, opS(ithCell,ithTimeBin,ithTraj,ithTrial)
%   - modelParam: the model parameter structure
%
% Outputs:
%   - linfields: linear place field cell array matching Jadhav lab
%       structure
%
% Optional:
%   - day, epoch, tetrode, plotAllFigs
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'opS',	@islogical)
addRequired(inputObj, 'modelParam',	@isstruct)
addParameter(inputObj, 'plotAllFigs', false, @islogical)
addParameter(inputObj, 'day', 1, @isnumeric)
addParameter(inputObj, 'epoch', 1, @isnumeric)
addParameter(inputObj, 'tetrode', 1, @isnumeric)
parse(inputObj, opS, modelParam, varargin{:});
p = inputObj.Results;


%% Set up

xPos = modelParam.xPos;
winSize = modelParam.linFieldGaussSD/modelParam.spatialBin*modelParam.winScale; % window is 5x the STDev
linTrajBins = numel(modelParam.gridxvals); % W-track: ~40 for arms, ~20 for half of top width, 2 cm bins


%% Calculate linfields data structure to match Jadhav lab data

% Calculate occupancy for each trajectory
accumOcc = zeros(linTrajBins, size(opS, 3));        % sum of occupancy timesteps along the 4 linear trajectories
for i = 1:numel(modelParam.t_PF)
    [~, linI] = min(( abs( [modelParam.gridxvals] - [xPos(i)] ) ));
    accumOcc(linI, :) = accumOcc(linI,:) + 1;
end
nTrials = size(opS, 4); % Don't assume nTrials==modelParam.nTrials_FP, in case a subset of trial spikes are passed in
accumOcc = accumOcc.*modelParam.dt*nTrials;

% Calcualte and store each cells place field data
for ithCell = 1:modelParam.n % Increment through cells

    accumTrajSpikes = zeros(linTrajBins, size(opS, 3)); % Initialize: sum of ithCell's spikes along the 4 linear trajectories
    for ithTraj = 1:size(opS, 3)
        spikeTimes = find( squeeze(sum(opS(ithCell,:,ithTraj,:), 4)) )*modelParam.dt;

        for i = 1:numel(spikeTimes)
            curX = xPos(round(spikeTimes(i)/modelParam.dt));        % xPos at spike time
            %[xL, linXind] = min(abs(modelParam.gridxvals - curX));  % Linear xPos index
            [~, linI] =     min(abs(modelParam.gridxvals - curX));    %
            accumTrajSpikes(linI, ithTraj) = accumTrajSpikes(linI,ithTraj)+1;
        end

        % Calculate linfields values:
        occupancy = accumOcc(:,ithTraj);         % Occupancy
        smoothedOccupancy = smoothdata(accumOcc(:,ithTraj),'gaussian',winSize); % smoothed occupancy
        spikeCount = accumTrajSpikes(:,ithTraj); % Spike count
        smoothedSpikeCount = smoothdata(accumTrajSpikes(:,ithTraj),'gaussian',winSize); % Smoothed spike count

        % Store linfields values
        linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,1) = [0:modelParam.spatialBin:(linTrajBins*modelParam.spatialBin-modelParam.spatialBin)]'*100; % (cm) linear distance from center well
        linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,2) = occupancy;
        linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,3) = spikeCount;
        linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,4) = spikeCount./occupancy; % occupancy-normalized firing rate
        linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,5) = smoothedSpikeCount./smoothedOccupancy; % smoothed occupancy-normalized firing rate
        linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,6) = smoothedOccupancy;
        linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,7) = smoothedSpikeCount;
    end
end


%% Plot some example figures, from the last cell

if p.plotAllFigs
    figure; plot(accumOcc); legend(); xlabel('Position bins'); ylabel('Occupancy')
    figure; plot(smoothdata(accumOcc,'gaussian',winSize)); legend(); xlabel('Position bins'); ylabel('Smoothed Occupancy')
end

if p.plotAllFigs
    figure; plot(accumTrajSpikes); legend(); xlabel('Position bins'); ylabel('Accumulated spike locations')
    figure; plot(smoothdata(accumTrajSpikes,'gaussian',winSize)); legend({'tr1', 'tr2', 'tr3', 'tr4'}); xlabel('Position bins'); ylabel('Smoothed spike count')
end

if p.plotAllFigs
    figure; plot(linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,1), linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,2))
    figure; plot(linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,1), linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,3))
    figure; plot(linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,1), linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,4))
    figure; plot(linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,1), linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,5))
    figure; plot(linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,1), linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,6))
    figure; plot(linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,1), linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,7))
end

end
