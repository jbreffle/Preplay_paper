function [figHandle, pValMat] = makePvalMatrix(eventRs, shuffleRs, eventJD, shuffleJD, varargin)
% [figHandle, pValMat] = utils.makePvalMatrix();
%
% Generates the preplay p-value matrix, as in Farooq et al., 2019.
%
% Input:
%   eventRs, column vector of actual event weighted correlations
%   shuffleRs, column cell array of shuffled event weighted correlations
%   eventJD, column vector of actual event maximum jump distances
%   shuffleJD, column cell array of shuffled event maximum jump distances
%
% Ouput:
%   figHandle, figure handle for the main p-value matrix
%   pValMat, the p-value matrix plotted in the figHandle figure
%


%% Parse inputs and set up for analysis
inputObj = inputParser;
addRequired(inputObj, 'eventRs',	@ismatrix)
addRequired(inputObj, 'shuffleRs',	@iscell)
addRequired(inputObj, 'eventJD',	@ismatrix)
addRequired(inputObj, 'shuffleJD',	@iscell)
addParameter(inputObj, 'plotExtra', false,	@islogical)
addParameter(inputObj, 'nanThreshIs0', true,	@islogical)
parse(inputObj, eventRs, shuffleRs, eventJD, shuffleJD, varargin{:});
p = inputObj.Results;

assert(iscolumn(eventRs))
assert(isequal(size(eventRs), size(shuffleRs), size(eventJD), size(shuffleJD)))

rvalThresh_vec = 0.0:0.1:0.9;
jumpThres_vec =  0.1:0.1:1.0;

pValMat = zeros(numel(jumpThres_vec), numel(rvalThresh_vec));
frac_shuff_events = zeros(numel(jumpThres_vec), numel(rvalThresh_vec));
frac_actu_events = zeros(numel(jumpThres_vec), numel(rvalThresh_vec));

nEvents = size(shuffleRs, 1);
nShuffles = size(shuffleRs{1}, 1);

if p.nanThreshIs0
    nanThresh = 0;
else
    nanThresh = (1/nEvents);
end

%% Calculate p-value matrix

for ithJD = 1:numel(jumpThres_vec)
    for ithR = 1:numel(rvalThresh_vec)

        fracEventPass = mean( (abs(eventRs)>rvalThresh_vec(ithR)) & (eventJD<jumpThres_vec(ithJD)) );

        fracShuffPass = zeros(1, nShuffles);
        for ithShuf = 1:nShuffles
            nShuff_jump = cellfun(@(x)     x(ithShuf) < jumpThres_vec(ithJD),  shuffleJD);
            nShuff_rval = cellfun(@(x) abs(x(ithShuf))> rvalThresh_vec(ithR), shuffleRs);
            fracShuffPass(ithShuf) = mean( nShuff_jump & nShuff_rval );
        end

        pValMat(ithJD, ithR) = 1 - mean(fracEventPass>fracShuffPass);
        if (mean(fracShuffPass)<=nanThresh) && (fracEventPass==0)
            pValMat(ithJD, ithR) = nan;
        end

        % Fraction of events meeting thresholds, for p.plotExtra
        frac_shuff_events(ithJD, ithR) = mean(fracShuffPass);
        frac_actu_events(ithJD, ithR) = fracEventPass;

        if jumpThres_vec(ithJD)==0.5 && rvalThresh_vec(ithR)==0.5

            plotExampleShuffle = false;
            if plotExampleShuffle
                exampleIthShuffle = 1;
                exampleShuffleJD = cellfun(@(x) x(exampleIthShuffle), shuffleJD);
                exampleShuffleR = cellfun(@(x) x(exampleIthShuffle), shuffleRs);

                figure; hold on; ecdf(eventJD); ecdf(exampleShuffleJD);
                legend('Actual', 'Shuffle', 'location', 'best')
                title('Comparison of one example shuffle CDF')
                xlabel('JD'); ylabel('CDF');

                figure; hold on; ecdf(eventJD); ecdf(exampleShuffleR);
                legend('Actual', 'Shuffle', 'location', 'best')
                title('Comparison of one example shuffle CDF')
                xlabel('|r|'); ylabel('CDF');
            end

            if p.plotExtra
                figure; hold on; xline(fracEventPass, 'k', 'LineWidth', 1); histogram(fracShuffPass);
                leg1 = legend('Simulation', 'Shuffle', 'location', 'best', 'box', 'off');
                leg1.ItemTokenSize = [10, 4];
                xlabel('Fraction of events'); ylabel({'Count', '(shuffle datasets)'});
            end

        end

    end
end

%% Plot results

% there is a strange bug where using exportgraphics() to save this figure
% as a pdf causes white squares to be saved as black. With this colormap
% white is exactly 0.05.
pValMat(abs(pValMat-0.05)<0.002) = 0.052;

% Plot with better colormap
figHandle = figure;
imagesc(jumpThres_vec, rvalThresh_vec, log10(pValMat'), 'AlphaData', ~isnan(pValMat'))
set(gca,'color',0*[1 1 1])
xlabel('<|Max Jump Distance|')
ylabel('>|Correlation|')
N = 256; n = N/2;
cm = NaN(N,3);
cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)'];
cm(:,3) = [linspace(0,1,n)';ones(N-n,1)];
set(gca,'clim',[log10(.05)*2 0])
set(gcf,'colormap',cm)
cb = colorbar('Direction','reverse','Ticks',[log10(.005),log10(.05),log10(.5)],'TickLabels',[.005,.05,.5]);
cb.Label.String = 'P-value';

if p.plotExtra
    % Plot with linear colormap
    figure; imagesc(jumpThres_vec, rvalThresh_vec, pValMat', 'AlphaData', ~isnan(pValMat')); colorbar
    set(gca,'color',0*[1 1 1])
    xlabel('<|Max Jump Distance|')
    ylabel('>|Correlation|')
    caxis([0, 0.1])

    % Plot mean num shuffled events passsing
    figure; imagesc(jumpThres_vec, rvalThresh_vec, frac_shuff_events', 'AlphaData', ~isnan(frac_shuff_events')); colorbar
    caxis([0, 1]); title('Mean frac shuff passing')
    xlabel('<|Max Jump Distance|')
    ylabel('>|Correlation|')

    % Plot mean num shuffled events passsing
    figure; imagesc(jumpThres_vec, rvalThresh_vec, frac_actu_events', 'AlphaData', ~isnan(frac_actu_events')); colorbar
    caxis([0, 1]); title('Mean frac actual events passing')
    xlabel('<|Max Jump Distance|')
    ylabel('>|Correlation|')
end

end
