function G_in = create_G_in(modelParam, network, ithEnv, isPFSim, varargin)
% G_in = create_G_in(modelParam, network, ithEnv, isPFSim)
%
% Creates a Poisson input conductance for either PF or preplay simulation.
%
% Inputs:
%   - modelParam: Structure containing model parameters
%   - network: Structure containing network data
%   - ithEnv: Index of the environment
%   - isPFSim: Flag indicating if it's a place field simulation (otherwise
%       a preplay simulation)
%
% Optional Input Parameters (parameter-value pairs):
%   - initialize (default: true): Flag indicating if G_in
%       should be initialized to steadystate values
%   - runTets (default: true): Flag indicating assertions
%       testing the code should be run (~15% of runtime due to prctile)
%
% Output:
%   - G_in: Matrix of Poisson input conductance values
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'network',    @isstruct)
addRequired(inputObj, 'ithEnv',     @isnumeric)
addRequired(inputObj, 'isPFSim',    @islogical)
addParameter(inputObj, 'initialize', true,	@islogical);
addParameter(inputObj, 'runTests', true,    @islogical);
addParameter(inputObj, 'seed', -1,    @isnumeric);
parse(inputObj, modelParam, network, ithEnv, isPFSim, varargin{:});
p = inputObj.Results;

% Only set the rng if a seed was passed in
if p.seed~=-1
    rng(p.seed, 'twister');
end

if isfield(modelParam, 'spatialInputType')
    inputType = modelParam.spatialInputType;
    if isequal(modelParam.spatialInputType, 'stepped')
        inputNSteps = modelParam.inputNSteps; % 5
    end
else % linear or stepped
    inputType = 'linear';
end

%% Set up
if isPFSim
    simDuration = modelParam.t_max_PF;
else
    simDuration = modelParam.t_max_preplay;
end
nTimeSteps = (simDuration/modelParam.dt)+1;

% Binary indices of excitatory and inhibitory neurons
Einds = ismember(1:modelParam.n, network.E_indices)';
Iinds = ismember(1:modelParam.n, network.I_indices)';

% Initialize G_in
G_in = zeros(modelParam.n, nTimeSteps);
if p.initialize
    G_in(:,1) = 1/2 * modelParam.Win_mean * 2 * modelParam.rG * modelParam.tau_syn_E + ...
        sqrt(1/2*modelParam.tau_syn_E*modelParam.Win_mean.^2*2*modelParam.rG).*randn(modelParam.n, 1) ;
end

% Pre-compute:
expDecayE = exp(-modelParam.dt/modelParam.tau_syn_E);
dtBaseRate = modelParam.dt* modelParam.rG;
if isPFSim % If place field simulation
    contextScalingE = Einds.*network.contextInput(:,ithEnv).* modelParam.PFcontextFrac;
    contextScalingI = Iinds.*network.contextInput(:,ithEnv).* modelParam.IcueScale_PF;
    spatialScaling1 = network.spatialInput{1}(:,ithEnv).*(1-modelParam.PFcontextFrac);
    spatialScaling2 = network.spatialInput{2}(:,ithEnv).*(1-modelParam.PFcontextFrac);
    if isequal(inputType, 'linearTernary')
        spatialScaling3 = network.spatialInput{3}(:,ithEnv).*(1-modelParam.PFcontextFrac);
    end
else % If preplay simulation
    if isfield(modelParam, 'useSleepContext') && modelParam.useSleepContext
        sleepContext = network.sleepContext;
    else
        sleepContext = network.contextInput(:,ithEnv);
    end
    contextScalingE = Einds.*sleepContext.* 1;
    contextScalingI = Iinds.*sleepContext.* modelParam.IcueScale ;
end

% Define functions that return the spatially modulated spike inputs
% Select one of several to use as an optional input to create_G_in
% "i" is the simulation time dt index

% Linearly scales across space, to peak at either end of track
linearInput = @(i) ...
    spatialScaling1.* (rand(modelParam.n, 1) < (dtBaseRate*(i/nTimeSteps))) + ...
    spatialScaling2.* (rand(modelParam.n, 1) < (dtBaseRate*((nTimeSteps-i)/nTimeSteps)));
% Three linearly varying inputs, third peaks at center
linearTernaryInput = @(i) ...
    spatialScaling1.* (rand(modelParam.n, 1) < (dtBaseRate* max(0, ((2*i)/(nTimeSteps)-1))              )) + ...
    spatialScaling2.* (rand(modelParam.n, 1) < (dtBaseRate* max(0, (2*(nTimeSteps-i)/nTimeSteps-1))     )) + ...
    spatialScaling3.* (rand(modelParam.n, 1) < (dtBaseRate* (1-(1/nTimeSteps*2)*abs(i-(nTimeSteps/2)))  ));
% Stepped input
steppedInput = @(i) spatialScaling1.* (rand(modelParam.n, 1)<(dtBaseRate*round(inputNSteps*i/nTimeSteps)/inputNSteps)) + ...
    spatialScaling2.* (rand(modelParam.n, 1)<(dtBaseRate*round(inputNSteps*(nTimeSteps-i)/nTimeSteps)/inputNSteps));

switch inputType
    case 'linear'
        spatialInputStep = linearInput;
    case 'linearTernary'
        spatialInputStep = linearTernaryInput;
    case 'stepped'
        spatialInputStep = steppedInput;
    otherwise
        error("create_G_in: unknown inputType option")
end

if ~isPFSim
    spatialInputStep = @(i) 0;
end

contextInputStep = @() contextScalingE .* (rand(modelParam.n, 1)<dtBaseRate) + ...
    contextScalingI .* (rand(modelParam.n, 1)<dtBaseRate);


%% Simulate G_in as a Poisson process

for tInd = 2:nTimeSteps

    G_in_step = G_in(:,tInd-1);
    % Exponential decay from last time step
    G_in_step = G_in_step*expDecayE;
    % Add spikes from each input source:
    G_in_step = G_in_step + spatialInputStep(tInd) + contextInputStep();
    G_in(:,tInd) = G_in_step;

end


%% Validate initialization
if p.initialize && p.runTests % prctile is ~ 15% of function call run-time
    initVals = G_in(:,1);
    simVals = G_in(:,2:end);
    P = prctile(simVals,[5, 95],2);
    fracInitsCentered = mean(initVals>P(:,1) & initVals<P(:,2));

    % At least 10% of inputs' initial values should fall within the simulated 90% CI
    % 10% is a very liberal lower bound, just to make sure it didn't break completely
    assert(fracInitsCentered>0.1);
    % Plots to check if something went wrong
    %{
    figure; plot([1:nTimeSteps], G_in(1:5,:));
    figure; hold on;
    histogram(simVals,  'facealpha',.5,'edgecolor','none', 'Normalization', 'pdf');
    histogram(initVals, 'facealpha',.5,'edgecolor','none', 'Normalization', 'pdf');
    legend({'Sim', 'Init'}, 'Location', 'Best')
    %}
end


end