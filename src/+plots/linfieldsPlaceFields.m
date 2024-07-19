function [] = linfieldsPlaceFields(linfields, network, modelParam, varargin)
% plots.linfieldsPlaceFields(linfields, network, modelParam)
%
% Plots the place fields in the linfields structure
%
% Inputs:
%   - linfields: linear place field cell array matching Jadhav lab
%       structure
%   - PFpeaksSequence: PFpeaksSequence{ithTraj} is a vector of the sequence
%       of PF peaks of the E-cells
%   - network: the network structure
%   - modelParam: the model parameter structure
%
% Outputs:
%   - No output. Generates plots.
%
% Optional:
%   - day, epoch, tetrode, plotNormRates, minPeakRate
%
% TODO: plot histogram of I cells PFs and peak rates?
% TODO: Add option to select which plots to generate
%   e.g. {'individual', 'histograms', 'cross-sorted', 'peak-corr'}
%

%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'linfields',          @iscell)
addRequired(inputObj, 'network',            @isstruct)
addRequired(inputObj, 'modelParam',         @isstruct)
addParameter(inputObj, 'day',       1, @isnumeric)
addParameter(inputObj, 'epoch',     1, @isnumeric)
addParameter(inputObj, 'tetrode',   1, @isnumeric)
addParameter(inputObj, 'plotNormRates',  true, @islogical);
addParameter(inputObj, 'minPeakRate',    0,    @isnumeric);
parse(inputObj, linfields, network, modelParam, varargin{:});
p = inputObj.Results;


[allEnvPFs, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network);


%% Plot place fields by self-peak order
for ithEnv = 1:modelParam.nEnvironments
    posBins = linfields{p.day}{p.epoch}{p.tetrode}{1}{ithEnv}(:,1);
    PFmat = allEnvPFs{ithEnv};
    PFmat_E = PFmat(network.E_indices,:);

    row_all_zeros = find(all( PFmat_E==0, 2)) ;
    row_n_all_zeros = find(~all( PFmat_E==0, 2)) ;

    [~,peakRateLocation] = max(squeeze(PFmat_E(row_n_all_zeros,:)), [], 2);
    [~,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    curr_PFpeaksSequence = [row_n_all_zeros(sortedCellIndsbyPeakRateLocation); row_all_zeros];

    peakRates = max(PFmat_E(curr_PFpeaksSequence,:), [], 2);

    if p.plotNormRates
        rateDenom1 = peakRates;
        caxmax = 1;
    else
        rateDenom1 = 1;
        caxmax = max(PFmat_E, [], 'all');
    end

    % For each trajectory, plot the PFs and the distribution of peaks
    figure; tiledlayout(1,2); sgtitle(['Env ID ', num2str(ithEnv)], fontweight='normal', fontsize=12)
    nexttile
    imagesc(posBins, 1:modelParam.n_E, PFmat_E(curr_PFpeaksSequence,:)./rateDenom1 );
    colorbar; caxis([0, caxmax])
    xlabel('Position (cm)'); ylabel('Cell (sorted)');
    nexttile
    histogram(peakRates)
    xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');
end


%% Plot place fields by all different peak sequences
if modelParam.nEnvironments>1
    figure
    for ithEnv = 1:modelParam.nEnvironments
        posBins = linfields{p.day}{p.epoch}{p.tetrode}{1}{ithEnv}(:,1);
        tmp_PFpeaksSequence = PFpeaksSequence{ithEnv};

        for jthEnv = 1:modelParam.nEnvironments

            PFmat = allEnvPFs{jthEnv};
            PFmat_E = PFmat(network.E_indices,:);
            if  mod(jthEnv,2) ~= mod(ithEnv,2)% mod(ithEnvSequence,2)==0 && ~[ithEnvData==ithEnvSequence]
                PFmat_E = fliplr(PFmat_E);
            end

            if p.plotNormRates
                rateDenom1 = max(PFmat_E(tmp_PFpeaksSequence,:), [], 2);
                caxmax = 1;
            else
                rateDenom1 = 1;
                caxmax = max(PFmat_E, [], 'all');
            end

            % Plot subplot's data
            ithPlot = sub2ind( [modelParam.nEnvironments, modelParam.nEnvironments], ithEnv, jthEnv);
            subplot(modelParam.nEnvironments, modelParam.nEnvironments, ithPlot);
            imagesc(posBins, 1:modelParam.n_E, PFmat_E(tmp_PFpeaksSequence,:)./rateDenom1 );

            % Put correlation in title
            pixelCorrs = diag(corr(PFmat_E, allEnvPFs{ithEnv}(network.E_indices,:) ));
            % assert(numel(pixelCorrs)==numel(posBins))
            mapCorr = nanmean(pixelCorrs);

            title(['r=', num2str(mapCorr, '%0.2f')], fontweight='normal')

            if jthEnv==ithEnv
                nthEnv = ceil(jthEnv/2);
                if mod(jthEnv,2)==0; trajDir='right';
                else; trajDir='left'; end
                title(['Env', num2str(nthEnv), ' ', trajDir], fontweight='bold')
            end

            if ithEnv==1 && jthEnv==modelParam.nEnvironments
                xlabel('Position (cm)'); ylabel('Cell (column sort)');
            end
        end
    end
end


%% Plot peak sequence correlations
if modelParam.nEnvironments>1
    figure
    for ithEnv = 1:modelParam.nEnvironments
        PFmat = allEnvPFs{ithEnv};
        PFmat_E = PFmat(network.E_indices,:);
        [peakRate, ~] = max(PFmat_E, [], 2);
        PFpeaksSequence1 = PFpeaksSequence{ithEnv};
        PFpeaksSequence1(peakRate<p.minPeakRate) = nan;

        for jthEnv = 1:modelParam.nEnvironments
            PFmat = allEnvPFs{jthEnv};
            PFmat_E = PFmat(network.E_indices,:);
            [peakRate, ~] = max(PFmat_E, [], 2);
            PFpeaksSequence2 = PFpeaksSequence{jthEnv};
            PFpeaksSequence2(peakRate<p.minPeakRate) = nan;

            % Plot subplot's data
            ithPlot = sub2ind( [modelParam.nEnvironments, modelParam.nEnvironments], ithEnv, jthEnv);
            subplot(modelParam.nEnvironments, modelParam.nEnvironments, ithPlot);
            scatter(PFpeaksSequence1, PFpeaksSequence2, 5);

            % Put correlation in title
            [RHO,~] = corr(PFpeaksSequence1, PFpeaksSequence2, 'rows','complete');
            title(['r=', num2str(RHO, '%0.2f')], fontweight='normal')

            if jthEnv==ithEnv
                nthEnv = ceil(jthEnv/2);
                if mod(jthEnv,2)==0; trajDir='right';
                else; trajDir='left'; end
                title(['Env', num2str(nthEnv), ' ', trajDir], fontweight='bold')
            end
            if ithEnv==1 && jthEnv==modelParam.nEnvironments
                xlabel('Peak rank'); ylabel('Peak rank');
            end
        end
    end
end


%% Finish
drawnow
myPlotSettings

end