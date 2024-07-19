function PFScore = calculate_linfieldsScore(linfields, modelParam, network, varargin)
% calculate_linfieldsScore
%
% Calculates the PF objective score, lower is better.
%
% Inputs:
%   - linfields: linear place field cell array matching Jadhav lab
%       structure
%   - modelParam: the model parameter structure
%   - network: the network structure
%
% Outputs:
%   - PFScore: the place field score for one environment/trajectory
%
% Optional:
%   - day, epoch, tetrode, trajectory
%

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'linfields',	@iscell)
addRequired(inputObj, 'modelParam',        @isstruct)
addRequired(inputObj, 'network',    @isstruct)
addParameter(inputObj, 'day', 1, @isnumeric)
addParameter(inputObj, 'epoch', 1, @isnumeric)
addParameter(inputObj, 'tetrode', 1, @isnumeric)
addParameter(inputObj, 'trajectory', 1, @isnumeric)
parse(inputObj, linfields, modelParam, network, varargin{:});
p = inputObj.Results;


% Calculates score on trajectory=1 if trajectory is not passed
if ~isempty(inputObj.UsingDefaults) && any(ismember(inputObj.UsingDefaults, 'trajectory'))
    if ~(numel(linfields{p.day}{p.epoch}{p.tetrode}{1})==1)
        warning(['linfields has multiple trajectories but the desired ', ...
            'trajectory to score was not specified. Scoring the first trajectory.'])
    end
end


%% Set up

% Unconstrained gaussian leads to high, persistant rates
% Constrain gaussian fit parameters to realistic place field properties
% 3 < peak firing rate < 30 Hz
% 0 < mean firing rate position < 100 cm
% 0 < standard deviation < 10 ~= firing above 25% of max covers ~25% of track

% gaus = @(x, a1, b1, c1,vo) a1*exp(-((x-b1)/c1).^2);
% figure; plot(posBins, gaus(posBins, 100, 50, 10) ); hold on; yline(25)
% a1 = amplitude = peak firing rate 3<peak<30
% b1 = mean = location of peak 0<mean<100 cm
% c1 = standard deviation = width of place field STD<10

% Defaults for TolFun and TolX are 1e-6
gaussFO = fitoptions('Method','NonlinearLeastSquares', ...
    'Lower',	modelParam.gaussFOLower, ...
    'Upper',	modelParam.gaussFOUpper, ...
    'TolFun',   1e-3, ...
    'TolX',     1e-3 ...
    ); % [peak amplitude, position of peak on track, standard deviation of peak]

posBins = linfields{p.day}{p.epoch}{p.tetrode}{1}{p.trajectory}(:,1);


%% Fit each PF to a constrained Gaussian

maxFR = nan(1, numel(network.all_indices));
allPFs = nan(numel(posBins), numel(network.all_indices));
gaussRsqrs = nan(1, numel(network.all_indices));
for ithCell = network.E_indices
    PF = linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{p.trajectory}(:,5);
    allPFs(:,ithCell)= PF;
    maxFR(ithCell) = max(PF);

    include=~isnan(PF);
    [~,gof,~] = fit(posBins(include), PF(include), 'gauss1', gaussFO); % can choose up to gauss8, for sum of 8 gaussians
    gaussRsqrs(ithCell) = gof.rsquare;
end


%% Calculate objective score, based on PF similarities to Gaussian

% Set all improper gaussRsqrs to minimum
gaussRsqrs(isnan(gaussRsqrs))= -1;
gaussRsqrs(isinf(gaussRsqrs))= -1;
gaussRsqrs(gaussRsqrs<-1) = -1;

% If a cell never fires, make it nan to exclude from analysis
maxFR(maxFR==0)=nan;

PFScore = 0;
for i = 1:numel(posBins)
    eCellInds = ismember(network.all_indices, network.E_indices);
    Inds = logical( ismember(network.all_indices, network.E_indices).* [allPFs(i,:)==maxFR] ); % get index of cells whose peak is at spatial bin i
    locMax = max(gaussRsqrs(logical(Inds))); % Determine best gaussian correlation of out of Inds cells

    if ~isempty(locMax)
        MeasureOfhighFiring = (1- [abs( maxFR(Inds)-modelParam.peakTarget ) / modelParam.peakTarget]); %index for proximity of a peak of 15 Hz
        PFScore = PFScore + -max( gaussRsqrs(Inds)+MeasureOfhighFiring ); % accumulate to the score the max gaussRsqr out of cells with peak at position i
    end
end
PFScore = PFScore/numel(posBins);

end