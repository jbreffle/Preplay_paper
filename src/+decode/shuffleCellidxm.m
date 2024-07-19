function shuffledCellidxm = shuffleCellidxm(cellidxm, modelParam, network)
% shuffledCellidxm = decode.shuffleCellidxm(cellidxm modelParam, network);
%
% Shuffles rows of the cellidxm matrix following different methods
%

switch modelParam.decodeCellShuffle

    case "clusterIndependent"
        % Shuffle independently of cluster membership
        hpnum = size(cellidxm, 1);
        shuffleInds = randperm(hpnum, hpnum);
        shuffledCellidxm = cellidxm(shuffleInds,:);

    case "withinCluster"
        % Shuffle only within clusters
        % Shuffle spikes by shuffling cell identity in cellidxm
        % 1) increment through cells and select a single cluster to
        % shuffle within for each
        validCellClusterMat = network.cluster_mat(:,network.E_indices(cellidxm(:,2)));
        shuffleClusterMembership = zeros(1, size(validCellClusterMat, 2));
        for ithValidCell = 1:numel(shuffleClusterMembership)
            candidates = find(validCellClusterMat(:,ithValidCell));
            shuffleClusterMembership(ithValidCell) = datasample(candidates, 1);
        end
        % 2) increment through clusters and shuffle the cells
        % based on shuffleClusterMembership
        shuffledCellidxm = nan(size(cellidxm));
        for ithCluster = 1:size(network.cluster_mat, 1)
            clusterCellInds = find(shuffleClusterMembership==ithCluster);
            shuffledClusterCellInds = clusterCellInds(randperm(numel(clusterCellInds), numel(clusterCellInds)));
            shuffledCellidxm(clusterCellInds,:) = cellidxm(shuffledClusterCellInds,:);
        end

    case "singleClusterCells"
        % Shuffle wtihin clusters only those cells in a single cluster
        % Same as withinCluster, but do not shuffle identity of cells
        % in multiple clusters
        % Shuffle spikes by shuffling cell identity in cellidxm
        % 1) increment through cells and select a single cluster to
        % shuffle within for each
        validCellClusterMat = network.cluster_mat(:,network.E_indices(cellidxm(:,2)));
        shuffleClusterMembership = nan(1, size(validCellClusterMat, 2));
        for ithValidCell = 1:numel(shuffleClusterMembership)
            cellCluster = find(validCellClusterMat(:,ithValidCell));
            if numel(cellCluster)>1
                cellCluster = nan;
            end
            shuffleClusterMembership(ithValidCell) = cellCluster;
        end
        % 2) increment through clusters and shuffle the cells
        % based on shuffleClusterMembership
        shuffledCellidxm = cellidxm;
        for ithCluster = 1:size(network.cluster_mat, 1)
            clusterCellInds = find(shuffleClusterMembership==ithCluster);
            shuffledClusterCellInds = clusterCellInds(randperm(numel(clusterCellInds), numel(clusterCellInds)));
            shuffledCellidxm(clusterCellInds,:) = cellidxm(shuffledClusterCellInds,:);
        end

    case "none"
        shuffledCellidxm = cellidxm;

    otherwise
        error("Unknown value for modelParam.decodeCellShuffle")
end