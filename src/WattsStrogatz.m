function h = WattsStrogatz(N, K, beta, varargin)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.
%
% Adapated from
% https://www.mathworks.com/help/matlab/math/build-watts-strogatz-small-world-graph-model.html


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'N',	@isnumeric)
addRequired(inputObj, 'K',        @isnumeric)
addRequired(inputObj, 'beta',    @isnumeric)
% addParameter(inputObj, 'runTests', true,    @islogical); % TODO: add optional plotPFs flag
parse(inputObj, N, K, beta, varargin{:});
p = inputObj.Results;


%% Main

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N
    switchEdge = rand(K, 1) < beta;

    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;

    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
end