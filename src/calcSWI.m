function [SWI, C, L] = calcSWI(W)
% CALCSWI calculates the Small-World Index of a connection matrix.
%
% Inputs:
%   - W: Connection matrix (weighted or binary)
%
% Output:
%   - SWI: Small-World Index of the matrix
%
% References:
% Neal et al., 2017
% https://www.cambridge.org/core/services/aop-cambridge-core/content/view/56ED3E109EE91CEA6FAE2D82FF20DB30/S2050124217000054a.pdf/how-small-is-it-comparing-indices-of-small-worldliness.pdf

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'W',	@islogical)
parse(inputObj, W);
p = inputObj.Results;


%% Main:

% Binarize the connection matrix
W = W > 0;

% Get the size of the matrix
n = size(W, 1);

% Calculate average degree
k = sum(W, 'all') / n;

% Calculate connection probability
p = sum(W, 'all') / (n^2 - n);

% Define Euler's gamma constant
egamma = 0.577215664901533; % equivalent to "double(eulergamma)" which is slow

% Calculate parameters for random networks (SWI right-hand side)
C_r = p;
L_r = (log(n) - egamma) / log(k) + 0.5;

% Calculate parameters for lattice networks (SWI left-hand side)
C_l = (3 * (k - 2)) / (4 * (k - 1));
L_l = (n) / (2 * k) + 0.5;

% Calculate pairwise distances
wDistances = distances(digraph(W));
wDistances(wDistances == 0) = nan; % Ignore self-distances with below nan-mean

% Calculate characteristic path length
L = nanmean(wDistances, 'all');

% Calculate clustering coefficient
C = calcFagioloCC(W);

% Calculate Small-World Index
SWI = (L - L_l) / (L_r - L_l) * (C - C_r) / (C_l - C_r);
end