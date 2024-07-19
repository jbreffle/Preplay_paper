%% Calculate SWI for randNet across nclusters and mnc
% swiGrid.m

addpath(fullfile("..", "..", "src"))
Config = utils.setConfig;

myPlotSettings()

% Select which panel to plot
panelToPlot = '8b'; % 1g, 8a, 8b, 8c

%% Main script parameters

switch panelToPlot
    case '1g'
        modelParam.conn_prob = 0.08; % Use {0.04, 0.08, 0.12} for Figure 4
        mncVec = 1:0.025:3.0; nClustersVec = 2:1:30; % Finer grid
    case '8a'
        modelParam.conn_prob = 0.04; % Use {0.04, 0.08, 0.12} for Figure 4
        mncVec = 1:0.25:3.0; nClustersVec = 5:5:30; % to match simulated grid
    case '8b'
        modelParam.conn_prob = 0.08; % Use {0.04, 0.08, 0.12} for Figure 4
        mncVec = 1:0.25:3.0; nClustersVec = 5:5:30; % to match simulated grid
    case '8c'
        modelParam.conn_prob = 0.012; % Use {0.04, 0.08, 0.12} for Figure 4
        mncVec = 1:0.25:3.0; nClustersVec = 5:5:30; % to match simulated grid
end

%modelParam.conn_prob = 0.08; % Use {0.04, 0.08, 0.12} for Figure 8
%mncVec = 1:0.25:3.0; nClustersVec = 5:5:30; % to match simulated grid
%mncVec = 1:0.025:3.0; nClustersVec = 2:1:30; % Finer grid

plotAllNets = false;
plotThreshScatter = false;
plotExtra = false;

plotThresholdLine = true;
plotBoundary = false;
MSPlots = true;

plotConnProb = false;
plotInterp = false;

nNets = 10;

modelParam.n = 500;


%%

% Options for plotting the boundary
SWIthres = 0.4;
SWIThreshType = fittype('a*exp(x) + c'); % Or without +c, fittype('exp1');%
startPoints = [3.5, -7];
%SWIThreshType = fittype('a*exp((x+b)/d) + c'); % Or without +c, fittype('exp1');%
% SWIThreshType = fittype('(a*x)^(b*2) + c'); % Or without +c, fittype('exp1');%
extendedMNCVec = (min(mncVec)*0.5):0.01:(max(mncVec)*1.5);
fullConnBoundary = @(mnc) (mnc.^2)/modelParam.conn_prob;

% Set up base network parameters
modelParam.p_E = 0.75; % Fraction of neurons that are excitatory
modelParam.p_I = 1-modelParam.p_E;
modelParam.n_I = round((1-modelParam.p_E)*modelParam.n);
modelParam.n_E = modelParam.n - modelParam.n_I;
modelParam.Win_mu = 1; modelParam.Win_sigma = 1; modelParam.IcueScale = 1;
modelParam.envIDs = 1;
modelParam.include_all = 3; % =3 means all neurons in at least one cluster
modelParam.global_inhib = true;


%%

tic
op_SWI_all = nan(numel(mncVec), numel(nClustersVec), nNets);
op_actualConnProb_all = nan(numel(mncVec), numel(nClustersVec), nNets);

rng('default')
for ithmnc = 1:numel(mncVec)
    modelParam.mnc = mncVec(ithmnc);

    for ithnClust = 1:numel(nClustersVec)
        modelParam.clusters = nClustersVec(ithnClust);

        if ~(fullConnBoundary(mncVec(ithmnc)) > nClustersVec(ithnClust))
            continue
        end
        for ithNet = 1:nNets
            %rng(ithNet)

            if modelParam.mnc>modelParam.clusters
                op_SWI(ithmnc, ithnClust) = nan;
                continue
            end

            modelParam.cluster_n = round((modelParam.mnc*modelParam.n) / modelParam.clusters);
            modelParam.npairs = modelParam.n*(modelParam.n-1);
            modelParam.nclusterpairs = modelParam.cluster_n*(modelParam.cluster_n - 1)*modelParam.clusters;
            modelParam.cluster_prob = min(modelParam.conn_prob*modelParam.npairs/modelParam.nclusterpairs,1);
            modelParam.envIDs=[];

            network = create_network(modelParam);
            W = network.conns(network.E_indices, network.E_indices);
            W = W>0;

            SWI_actual = calcSWI(W);

            op_SWI_all(ithmnc, ithnClust, ithNet) = SWI_actual;
            op_actualConnProb_all(ithmnc, ithnClust, ithNet) = mean(mean(W, 1), 2);

            if plotAllNets
                % figure; imagesc(W) % Plot connectivity matrix
                maxNClust = modelParam.clusters;
                x = W; p = squareform(pdist(x)); l = linkage(p); c = cluster(l, 'Maxclust',maxNClust); [~, I] = sort(c);
                figure; imagesc(x(I,I)); colorbar
                title(['Clustered W, nClust=', num2str(maxNClust)]); xlabel('Neuron'); ylabel('Neuron');
            end

        end % End net loop
    end     % End clusters loop
end         % End mnc loop

runTime = toc;
disp(['runtime: ', num2str(runTime)])


%% Calculate the median values of only valid networks

op_actualConnProb_all(isinf(op_actualConnProb_all))=nan;
op_SWI_all(isinf(op_SWI_all))=nan;
op_SWI = nanmedian(op_SWI_all, 3);
op_actualConnProb = nanmedian(op_actualConnProb_all, 3);

wrongConnectivity = ([(1./nClustersVec).*(modelParam.n_E.*mncVec')].*mncVec') < modelParam.n_E*modelParam.conn_prob;
badSWI = [isnan(op_SWI) | isinf(op_SWI)];
op2 = op_SWI;
op2(wrongConnectivity) = nan;
isNanVisbility = 1.0;
AlphaData =  ~isnan(op2) + [isnan(op2)*isNanVisbility];
AlphaData(badSWI) = 0;


%% Plot:

if MSPlots
    myPlotSettings(width=1.75, height=1.35)
end

figure; imagesc(mncVec, nClustersVec, op_SWI', 'AlphaData', AlphaData');
disp(['Max SWI: ', num2str(max(op_SWI, [], 'all'))])
title(['$p_c=', num2str(modelParam.conn_prob), '$'], 'FontWeight','Normal','interpreter','latex')
xlabel('Cluster participation'); ylabel('Clusters')
% halfPerctlEdge = 20; caxis([prctile(op, halfPerctlEdge, 'all'), prctile(op, 100-halfPerctlEdge, 'all')])
% caxis([0, 1.2]) % For grid matching Fig 3
cb = colorbar('XTick', [0,1]); set(gca,'YDir','normal')
cb.Label.String = 'SWI';
cb.Limits = [0, 1];
cb.Ticks = [0, 1];

% Plot the boundary where the network cannot match parameters.conn_prob
if plotBoundary
    hold on; plot(extendedMNCVec, fullConnBoundary(extendedMNCVec), 'k:', LineWidth=1.5)
end

% Detect threshold crossing, then plot fitted exponential curve
if plotThresholdLine
    x = op2>SWIthres;
    % x(isnan(op2)) = nan;
    [xvalinds, yvalinds] = find( diff(x, [], 2)==1 ) ;
    [fexp, fexp_G] = fit(mncVec(xvalinds)', nClustersVec(yvalinds)', SWIThreshType, 'StartPoint', startPoints);
    hold on; plot(extendedMNCVec, fexp(extendedMNCVec), 'r:', LineWidth=1.5)
    fexp
    if plotThreshScatter
        scatter(mncVec(xvalinds), nClustersVec(yvalinds), 'r', LineWidth=1.5)
    end
    % figure; scatter(mncVec(xvalinds)', nClustersVec(yvalinds)')
end

% Plot the actual mean E-E connection probability
if plotConnProb
    figure; imagesc(mncVec, nClustersVec, op_actualConnProb', 'AlphaData', AlphaData');
    colorbar; set(gca,'YDir','normal')
    %disp(['Max SWI: ', num2str(max(op_SWI, [], 'all'))])
    title(['Conn. prob: p_{c}=', num2str(modelParam.conn_prob)], 'FontWeight','Normal')
    xlabel('Cluster participation'); ylabel('Clusters')
    caxis([modelParam.conn_prob*0.5, modelParam.conn_prob*1])
end

if 0
    % Difference in SWI, if op16=op_SWI for p=0.16, and op08=op_SWI for p=0.08
    figure; imagesc(mncVec, nClustersVec, (op16-op08)', 'AlphaData', AlphaData');
    title('Difference in SWI')
    colorbar; set(gca,'YDir','normal')
    xlabel('Cluster participation'); ylabel('Clusters')
end

% Additional plots related to above
if plotExtra
    % How many neurons are in each cluster
    nNeuronsPerCluster = [(1./nClustersVec).*(modelParam.n_E.*mncVec')];
    figure; imagesc(mncVec, nClustersVec, nNeuronsPerCluster'); %, 'AlphaData', ~isnan(op'));
    colorbar; set(gca,'YDir','normal')
    title('nNeuronsPerCluster'); xlabel('Cluster participation'); ylabel('nClusters')

    % How many neurons can each neuron connect to
    nPostSynCand = nNeuronsPerCluster.*mncVec';
    figure; imagesc(mncVec, nClustersVec, nPostSynCand'); %, 'AlphaData', ~isnan(op'));
    colorbar; set(gca,'YDir','normal')
    title('nPostSynCand'); xlabel('Cluster participation'); ylabel('nClusters')

    % What fraction of possible connections are present
    fracPostConn =  (modelParam.n_E*modelParam.conn_prob) ./nPostSynCand; fracPostConn(fracPostConn>1) = 1;
    figure; imagesc(mncVec, nClustersVec, fracPostConn'); %, 'AlphaData', ~isnan(op'));
    colorbar; set(gca,'YDir','normal')
    title('fracPostConn'); xlabel('Cluster participation'); ylabel('nClusters')
end

if plotInterp
    % TODO: use interpolation for fitting the threshold

    xIntrpVals = min(mncVec):mean(diff(mncVec))/10:max(mncVec);
    yIntrpVals = min(nClustersVec):mean(diff(nClustersVec))/10:max(nClustersVec);
    [xq,yq] = meshgrid(xIntrpVals, yIntrpVals);
    [xActual,yActual] = meshgrid(mncVec, nClustersVec);
    vq = griddata(xActual,yActual,op_SWI',xq,yq);

    % Imagesc plot
    figure; imagesc(xIntrpVals,yIntrpVals,vq)
    colorbar; set(gca,'YDir','normal')
    x = double(vq>SWIthres);
    x(isnan(vq))=nan;
    [xvalinds, yvalinds] = find( diff(x', [], 2)==1 ) ;
    [fexp_interp, fexp_G] = fit(xIntrpVals(xvalinds)', yIntrpVals(yvalinds)', SWIThreshType, 'StartPoint', startPoints);
    hold on; plot(extendedMNCVec, fexp_interp(extendedMNCVec), 'r:', LineWidth=1.5)
    % figure; scatter(xIntrpVals(xvalinds)', yIntrpVals(yvalinds)')
    fexp_interp

end