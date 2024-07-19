function [] = randnetDecode(replaytrajectory, varargin)
% plots.randnetDecode.m
%
% Plots decode results from run_net_decode.m
%
% Inputs:
%   - replaytrajectory: Bayesian decoding results cell array
%       - replaytrajectory{p.day}{p.ep}.<fieldname>
%
% Outputs:
%   - No output. Generates plots.
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'replaytrajectory',	@iscell)
addParameter(inputObj, 'day', 1, @isnumeric);
addParameter(inputObj, 'ep', 2, @isnumeric);
addParameter(inputObj, 'ithTraj', 1, @isnumeric);
addParameter(inputObj, 'plotDecodeDist', true, @islogical);
addParameter(inputObj, 'plotPvalMat', true, @islogical);
parse(inputObj, replaytrajectory, varargin{:});
p = inputObj.Results;


%% Plot preplay decoding results

if p.plotDecodeDist
    pvals = replaytrajectory{p.day}{p.ep}.pvalue(:,p.ithTraj);
    figure; histogram(pvals, 10)
    xlabel('P-value'); ylabel('Event (count)')
    title(['EnvID: ', num2str(p.ithTraj)])

    allshuff_rvals = vertcat(replaytrajectory{p.day}{p.ep}.shuffle_rsquare{:});
    allshuff_rvals = allshuff_rvals(:,p.ithTraj); % Takling just forward traj

    rvals_preplay = replaytrajectory{p.day}{p.ep}.rsquare(:,p.ithTraj);
    rvals_shuffle = vertcat(allshuff_rvals{:,1});
    figure; hold on;
    ecdf(rvals_preplay)
    ecdf(rvals_shuffle)
    xlabel('Decode abs(r)'); ylabel('CDF')

    [~,P,~] = kstest2(rvals_preplay, rvals_shuffle);
    disp(['Decode pval: ', num2str(P)])
    legend({'Preplays', 'Shuffles'}, 'Location', 'Best')
    title(['EnvID: ', num2str(p.ithTraj)])
end


%% rval and max jump threshold pval matrix
if p.plotPvalMat
    allshuff_rvals = vertcat(replaytrajectory{p.day}{p.ep}.shuffle_rsquare{:});
    allshuff_rvals = allshuff_rvals(:,p.ithTraj); % Takling just forward traj

    allshuff_jumps = vertcat(replaytrajectory{p.day}{p.ep}.shuffle_maxJump{:});
    allshuff_jumps = allshuff_jumps(:,p.ithTraj); % Takling just forward traj

    rvals_preplay = replaytrajectory{p.day}{p.ep}.rsquare(:,p.ithTraj);
    jump_preplay = replaytrajectory{p.day}{p.ep}.maxJump(:,p.ithTraj);

    rvalThresh_vec = 0:0.1:1;
    jumpThres_vec = 0:0.1:1;

    op = zeros(numel(jumpThres_vec), numel(rvalThresh_vec));
    for ij = 1:numel(jumpThres_vec)

        for ir = 1:numel(rvalThresh_vec)
            nActPass = mean( [rvals_preplay>rvalThresh_vec(ir)] & [jump_preplay<jumpThres_vec(ij)]);
            nShuffPass = zeros(1, size(allshuff_jumps{1}, 2));
            for ithShuf = 1:size(allshuff_jumps{1}, 2)
                nShuff_jump = cellfun(@(x) [x(ithShuf)<jumpThres_vec(ij)], allshuff_jumps);
                nShuff_rval = cellfun(@(x) [x(ithShuf)>rvalThresh_vec(ir)], allshuff_rvals);
                nShuffPass(ithShuf) = mean( nShuff_jump & nShuff_rval );
            end

            op(ij, ir) = 1 - mean(nActPass>nShuffPass);
            if (sum(nShuffPass)==0) && (sum(nActPass)==0)
                op(ij, ir) = nan;
            end
        end
    end

    %{
    figure; imagesc(rvalThresh_vec, jumpThres_vec, op', 'AlphaData', ~isnan(op')); colorbar
    xlabel('max jump')
    ylabel('r^2')
    %}

    % Plot with log-scale and 0.05 centered colormap
    figure;
    imagesc(rvalThresh_vec, jumpThres_vec, log10(op'), 'AlphaData', ~isnan(op'))
    % set(gca,'YDir','normal')
    cb = colorbar(); %cb.Label.String = cbLabel2;
    xlabel('Max jump')
    ylabel('abs(r)')
    %title(analysisTitle)

    N = 256; n = N/2;
    cm = NaN(N,3);
    cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
    cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)'];
    cm(:,3) = [linspace(0,1,n)';ones(N-n,1)];

    set(gca,'clim',[log10(.05)*2 0])
    set(gcf,'colormap',cm)
    colorbar
    colorbar('Direction','reverse','Ticks',[log10(.005),log10(.05),log10(.5)],'TickLabels',[.005,.05,.5])

    hold on;
    [xnan, ynan] = find(isnan(op));
    scatter(rvalThresh_vec(xnan), jumpThres_vec(ynan), 300, 'x')
    title(['EnvID: ', num2str(p.ithTraj)])
end


end