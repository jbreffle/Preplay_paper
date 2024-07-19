%% Main function to generate tests
%
% These tests can be run from the command window with "runtests myTests"
% while in the test folder, or with "runtests("test/myTests.m")" from the
% project folder
%
% Capture the results by saving the output: e.g. testResults = runtests("test/myTests.m")
%
% See different types of verify<X>() functions:
%   https://www.mathworks.com/help/matlab/matlab_prog/types-of-qualifications.html
%

function tests = myTests
addpath(fullfile("..", "src"))
Config = utils.setConfig;
tests = functiontests(localfunctions);
end


%% Test Functions

function test_param_structs(testCase)
[modelParam, simParam] = initialize_structs();
actualResult = isstruct(modelParam)&isstruct(simParam);
verifyTrue(testCase, actualResult)
end

function test_create_network(testCase)
[modelParam, ~] = initialize_structs();
network = create_network(modelParam);
actualResult = isstruct(network);
verifyTrue(testCase, actualResult)
end

function test_create_G_in(testCase)
[modelParam, ~] = initialize_structs();
network = create_network(modelParam);
ithEnv = 1;
isPFSim = true;
G_in = create_G_in(modelParam, network, ithEnv, isPFSim);
actualResult = size(G_in);
expected_result = [modelParam.n, modelParam.t_max_PF/modelParam.dt+1];
verifyEqual(testCase, actualResult, expected_result)
end

function test_create_G_in2(testCase)
rng('default')
[modelParam, ~] = initialize_structs();
network = create_network(modelParam);
ithEnv = 1;
isPFSim = true;
G_in = create_G_in(modelParam, network, ithEnv, isPFSim);
actualResult = abs((G_in(1)-5.304317916047892e-09))<1e-14 & ...
    abs((G_in(end)-4.486395830507146e-09))<1e-14;
verifyTrue(testCase, actualResult)
end

function test_create_G_in3(testCase)
rng('default')
[modelParam, ~] = initialize_structs();
network = create_network(modelParam);
ithEnv = 1;
isPFSim = true;
G_in = create_G_in(modelParam, network, ithEnv, isPFSim);
initVals = G_in(:,1);
simVals = G_in(:,2:end);
P = prctile(simVals,[5, 95],2);
fracInitsCentered = mean(initVals>P(:,1) & initVals<P(:,2));
actualResult = fracInitsCentered>0.1;
verifyTrue(testCase, actualResult)
end

function test_create_G_in4(testCase)
rng('default')
[modelParam, ~] = initialize_structs();
network = create_network(modelParam);
ithEnv = 1;
isPFSim = true;
G_in = create_G_in(modelParam, network, ithEnv, isPFSim, initialize=false);
initVals = G_in(:,1);
simVals = G_in(:,2:end);
P = prctile(simVals,[5, 95],2);
fracInitsCentered = mean(initVals>P(:,1) & initVals<P(:,2));
actualResult = fracInitsCentered>0.1;
verifyFalse(testCase, actualResult)
end


%% Helper functions

function [modelParam, simParam] = initialize_structs()
simParam.nNets          = 1;    % Number of networks to simulate at each parameter point
simParam.nTrials_preplay= 1;    % Number of preplay trials to simulate per network
simParam.t_max_preplay  = 10;   % Trial duration (s) for each preplay trial
simParam.nTrials_PF     = 5;    % Number of PF trials to simulate (each run across the track)
simParam.t_max_PF       = 2;    % Duration (s) for each PF trial
simParam.envIDs         = [1];	% Environmnent IDs to simulate for PFs (default is [1:4] if this is empty or non-existent)
simParam.calcPFScore    = false;
simParam.dispFlag_decode = false;

paramSet = 'PFtuned1' ;	% 'default', '7_22', '10_18', 'ISregMin', 'PFtuned1', 'PFtuned2'
modelParam = initialize_parameters(paramSet, simParam);
end
