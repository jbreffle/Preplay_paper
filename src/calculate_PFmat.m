function [PFmat, PFpeaksSequence] = calculate_PFmat(linfields, modelParam, network, varargin)
% calculate_PFmat
%
% Calculate the matrix of place fields PFmat and the sequence of place
% field peaks PFpeaksSequence from the linfields cell array
%
% Inputs:
%   - linfields: linear place field cell array matching Jadhav lab
%       structure
%   - modelParam: the model parameter structure
%   - network: the network structure
%
% Outputs:
%   - PFpeaksSequence: PFpeaksSequence{ithTraj} is a vector of the sequence
%       of PF peaks of the E-cells
%   - PFmat: PFmat{ithTraj} its the place field matrix for all cells
%       (includes I-cells)
%
% Optional:
%   - day, epoch, tetrode, plotPFsFlag
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'linfields',	@iscell)
addRequired(inputObj, 'modelParam',	@isstruct)
addRequired(inputObj, 'network',	 @isstruct)
addParameter(inputObj, 'plotPFsFlag', false, @islogical)
addParameter(inputObj, 'day', 1, @isnumeric)
addParameter(inputObj, 'epoch', 1, @isnumeric)
addParameter(inputObj, 'tetrode', 1, @isnumeric)
parse(inputObj, linfields, modelParam, network, varargin{:});
p = inputObj.Results;


%% Calculate the PF matrix PFmat and the sequence of PF peaks PFpeaksSequence
nTraj = numel(linfields{p.day}{p.epoch}{p.tetrode}{1});
PFmat = cell(1, nTraj);
PFpeaksSequence = cell(1, nTraj);

for ithTraj = 1:nTraj

    % Calculate the PF matrix
    PFmat_temp = [];
    for ithCell = 1:modelParam.n
        PF = linfields{p.day}{p.epoch}{p.tetrode}{ithCell}{ithTraj}(:,5);
        PFmat_temp = [PFmat_temp;PF'];
    end
    PFmat{ithTraj} = PFmat_temp;

    % Calculate the sequence of peaks
    PFmat_E = PFmat_temp(network.E_indices,:);
    row_all_zeros = find(all( PFmat_E==0, 2)) ;
    row_not_all_zeros = find(~all( PFmat_E==0, 2)) ;
    [~,peakRateLocation] = max(squeeze(PFmat_E(row_not_all_zeros,:)), [], 2);
    [~,sortedCellIndsbyPeakRateLocation] = sort(peakRateLocation, 'descend');
    PFpeaksSequence{ithTraj} = [row_not_all_zeros(sortedCellIndsbyPeakRateLocation); row_all_zeros];

    % Optionally, plot the place fields
    if p.plotPFsFlag
        [~, peakRateLocation_all] = max(PFmat_E, [], 2);
        normRates = 1;
        if normRates
            rateDenom1 = max(PFmat_E(PFpeaksSequence{ithTraj},:), [], 2);
            caxmax = 1;
        else
            rateDenom1 = 1;
            caxmax = max(PFmat_E, [], 'all');
        end
        % Place fields
        figure; imagesc( PFmat_E(PFpeaksSequence{ithTraj},:)./rateDenom1 ); title(['EnvID: ', num2str(ithTraj)]); colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
        % Place field peak rates
        figure; histogram(rateDenom1);  title(['EnvID: ', num2str(ithTraj)]);
        xlabel('Peak PF rate (Hz)'); ylabel('E cells (count)');
        % Place fields centered on peaks
        padPFmat = [PFmat_E, zeros( [size(PFmat_E, 1), size(PFmat_E, 2) ] )];
        shiftpadPFmat = cell2mat(arrayfun(@(i){ [circshift(padPFmat(i,:), [50-peakRateLocation_all(i)] )]' }, 1:numel(peakRateLocation_all)))' ;% output matrix
        figure; imagesc( shiftpadPFmat(PFpeaksSequence{ithTraj},:)./rateDenom1 ); title('PFs aligned to peak'); colorbar; caxis([0, caxmax])
        xlabel('Position (2 cm bin)'); ylabel('Cell (sorted)');
    end

end

end
