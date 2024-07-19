function network = create_network(modelParam, varargin)
% network = create_network(modelParam)
%
% Generates the network structure.
%
% Inputs:
%   - modelParam: Structure containing model parameters
%
% Optional:
%   - seed: random number generator seed
%
% Output:
%   - network: The network structure for a given net


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'parameters',	@isstruct)
addParameter(inputObj, 'seed', -1,    @isnumeric);
parse(inputObj, modelParam, varargin{:});
p = inputObj.Results;

% Only set the rng if a seed was passed in
if isempty(inputObj.UsingDefaults) || all(~ismember(inputObj.UsingDefaults, 'seed'))
    rng(p.seed, 'twister');
end

network = struct; % Initialize output


%% Create recurrent network

% Assign E/I cell identity
all_indices = [1:modelParam.n];
%I_indices = sort(datasample(all_indices,modelParam.n_I,'Replace',false)); % TODO: switch to this line
I_indices = datasample(all_indices,modelParam.n_I,'Replace',false);
E_indices = find(~ismember(all_indices,I_indices));

% Assign cluster membership
cluster_mat = zeros(modelParam.clusters,modelParam.n);
for i = 1:modelParam.clusters
    ord = randperm(modelParam.n,modelParam.cluster_n); % Random connectivity
    cluster_mat(i,ord) = 1;                            % Note which neurons are in a cluster
end
clear ord i
%cluster_mat = rand(modelParam.clusters,modelParam.n)<modelParam.cluster_memb_prob; % TODO: replace above with just this line

% TODO: remove options 1 and 2, only allow 3
if modelParam.include_all == 1
    %Add back in those neurons that are not placed in a cluster, by
    %removing a place from another neuron with a high presence - this
    %section can be removed if you'd like some neurons to be unconnected
    ind_non = find(sum(cluster_mat) == 0);
    ind_high = find(sum(cluster_mat) > 2);
    for i = ind_non
        clust_place = randi(modelParam.clusters);
        ind_inclust = find(cluster_mat(clust_place,:));
        try %#ok<TRYNC>
            val_steal = datasample(intersect(ind_inclust,ind_high),1);
            cluster_mat(clust_place,i) = 1;
            cluster_mat(clust_place,val_steal) = 0;
            ind_high = find(sum(cluster_mat) > 2);
        end
    end
elseif modelParam.include_all==2
    % Same as include_all==1, but ind_high threshold is changed to 1
    % Note: if mnc=1, then by chance there may still be some neurons
    % in 2 clusters (if by chance mean mnc of cluster mat > 1)
    ind_non = find(sum(cluster_mat) == 0);
    ind_high = find(sum(cluster_mat) > 1);
    for i = ind_non
        clust_place = randi(modelParam.clusters);
        ind_inclust = find(cluster_mat(clust_place,:));
        try %#ok<TRYNC>
            val_steal = datasample(intersect(ind_inclust,ind_high),1);
            cluster_mat(clust_place,i) = 1;
            cluster_mat(clust_place,val_steal) = 0;
            ind_high = find(sum(cluster_mat) > 1);
        end
    end
elseif modelParam.include_all==3
    % All neurons assigned to one cluster, then mnc-1 is the new mnc
    % for assigning membership to additional clusters.
    % Overwrites above created cluster_mat.
    % Assumes mnc>=1.
    firstClusts=zeros(modelParam.clusters, modelParam.n);
    initialClusterInds =  randi(modelParam.clusters, 1, modelParam.n);
    initialClusterInds = sub2ind(size(firstClusts), initialClusterInds, 1:modelParam.n);
    firstClusts(initialClusterInds) = 1;
    pMem = (modelParam.mnc-1)/(modelParam.clusters-1); % probability a given neuron is a member of a given (additional) cluster
    randClusts = [rand(modelParam.clusters, modelParam.n)<pMem];
    cluster_mat = firstClusts|randClusts; % Combine first and additional cluster membership matrices
    % mean(sum(cluster_mat,1))
end

% Create connections
conns = zeros(modelParam.n);
for i = 1:modelParam.clusters
    ord = find(cluster_mat(i,:));
    ord_len = length(ord);
    conns(ord,ord) = conns(ord,ord) + (rand(ord_len,ord_len) <= modelParam.cluster_prob);
end

% Remove self-connectivity
for i = 1:modelParam.n
    conns(i,i) = 0;
end

% If global_inhib, overwrite all inhibitory outputs and inputs
if modelParam.global_inhib
    conns(I_indices,:) = (rand([modelParam.n*(1-modelParam.p_E), modelParam.n]) < modelParam.p_I);
    conns(:,I_indices) = (rand(modelParam.n, [modelParam.n*(1-modelParam.p_E)]) < modelParam.p_I);
end


%% Sleep context inputs, for simulating preplay
% Identical distribution as contextInput, but unique draw of values for preplay simulation

% Note: creating newtork.sleepContext will affect rng calls below it
if isfield(modelParam, 'useSleepContext') && modelParam.useSleepContext
    network.sleepContext = lognrnd(modelParam.Win_mu, modelParam.Win_sigma/4, modelParam.n, 1);
end

%% Create environment inputs
for ithEnv = modelParam.envIDs

    if mod(ithEnv,2) == 0  % If envID is even,
        % Use same context input for reverse traversal
        network.contextInput(:,ithEnv) = network.contextInput(:,ithEnv-1);

        network.clusterSequence(:,ithEnv) = flipud(network.clusterSequence(:,ithEnv-1));

        % Flip inputs1 and 2
        network.spatialInput{1}(:,ithEnv) = network.spatialInput{2}(:,ithEnv-1);
        network.spatialInput{2}(:,ithEnv) = network.spatialInput{1}(:,ithEnv-1);

    else % If envID is odd, create new input strengths for the new environment

        % Input strengths
        input1 = lognrnd(modelParam.Win_mu, modelParam.Win_sigma, modelParam.n, 1); % location cue 1 strengths
        input2 = lognrnd(modelParam.Win_mu, modelParam.Win_sigma, modelParam.n, 1); % location cue 2 strengths
        contextInput = lognrnd(modelParam.Win_mu, modelParam.Win_sigma/4, modelParam.n, 1); % context cue strength
        % contextInput = parameters.Win_mean+(sqrt(parameters.Win_var)*randn(parameters.n, 1)); % context cue strength
        %contextInput = lognrnd(log(parameters.Win_mu), parameters.Win_sigma, parameters.n, pfsim.nEnvironments);

        % if parameters.clusterCorrs==1, then make input1 and input2 correlated with cluster membership
        if [isfield(modelParam, 'clusterCorrs') && modelParam.clusterCorrs]
            clusterSequence = randperm(modelParam.clusters);  % 1:parameters.clusters;
            for i = 1:modelParam.n
                secondBias =  mean(cluster_mat(:,i).*(clusterSequence'-1)) / (1-1/modelParam.clusters) / sum(cluster_mat(:,i));
                secondBias = 0.5 + ((1-(2*secondBias))/modelParam.inputBiasSigma);
                temp1 = (input1(i)+input2(i)) * (1-secondBias);
                temp2 = (input1(i)+input2(i)) * (secondBias);
                input1(i) = temp1;
                input2(i) = temp2;
            end
        else
            clusterSequence = nan(1, modelParam.clusters);
        end

        % I cells don't receive spatially modulated input
        input1(I_indices) = 0;
        input2(I_indices) = 0;

        network.spatialInput{1}(:,ithEnv) = input1;
        network.spatialInput{2}(:,ithEnv) = input2;
        network.clusterSequence(:,ithEnv) = clusterSequence;
        network.contextInput(:,ithEnv) = contextInput;
    end
end

% Add a third input cue, peaked at the center of the track
if isfield(modelParam, 'spatialInputType') && isequal(modelParam.spatialInputType, 'linearTernary')
    % Create input strengths for the third cue
    if isfield(modelParam, 'clusterCorrs') && modelParam.clusterCorrs
        network.spatialInput{3} = (network.spatialInput{1} + network.spatialInput{2})./2;
    else
        network.spatialInput{3} = lognrnd(modelParam.Win_mu, modelParam.Win_sigma, modelParam.n, modelParam.nEnvironments);
        % I-cells don't receive spatially modulated input
        network.spatialInput{3}(I_indices,:) = 0;
        % Left/right trajectories of same env should be identical
        for ithEnv = 2:2:modelParam.nEnvironments
            network.spatialInput{3}(:,ithEnv) = network.spatialInput{3}(:,ithEnv-1);
        end
    end
end


%% Store fields in network struct
network.cluster_mat = cluster_mat;
network.conns = conns;
network.seed = p.seed;
network.n_EE = sum(network.conns(E_indices,E_indices),'all'); %number of E-E connections
network.n_EI = sum(network.conns(E_indices,I_indices),'all'); %number of E-I connections
network.n_II = sum(network.conns(I_indices,I_indices),'all'); %number of I-I connections
network.n_IE = sum(network.conns(I_indices,E_indices),'all'); %number of I-E connections
network.all_indices = all_indices;
network.I_indices = I_indices;
network.E_indices = E_indices;
network.validationRNG = rand(1); % When recreating a net from network.seed, use this value to verify exact re-creation

end