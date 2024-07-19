function fcc = calcFagioloCC(W)
% CALCFAGIOLOCC calculates the clustering coefficient of a directed binary graph.
% The function uses equations 6-8 from Fagiolo (2007) to compute the clustering coefficient.
%
% Inputs:
%   - W: Directed binary graph (adjacency matrix)
%        W must be a square matrix with binary values (0 or 1)
%
% Output:
%   - fcc: Clustering coefficient of the graph

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'W',	@islogical)
parse(inputObj, W);
p = inputObj.Results;


%% Main:
% Calculate total degree of each node (eq 6)
dtot = sum(W' + W, 2);

% Calculate the square of the adjacency matrix (eq 7)
dbi = sum(W.*W',2); % Faster equivalent to dbi = diag(W^2);

% Calculate node clustering coefficients (eq 8)
WW = W+W';
numerator = sum(WW.*(WW*WW))'; % Faster equivalent to numerator = diag((W + W')^3);

denominator = 2 .* (dtot .* (dtot - 1) - 2 .* dbi);
nodeCCs = numerator ./ denominator;

% Calculate the mean clustering coefficient
fcc = mean(nodeCCs);
end