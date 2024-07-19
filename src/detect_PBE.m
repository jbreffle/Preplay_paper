function [PBEstruct] = detect_PBE(spikes, modelParam, varargin)
% detect_PBE.m
%
% Detects population burst events. Assumes 'spikes' only includes cells
% that you want to include in burst detection (e.g. E-cells only).
%
% Inputs:
%   - opS: logical spike matrix, opS(ithCell,ithTimeBin,ithTraj,ithTrial)
%   - modelParam: the model parameter structure
%
% Output:
%   - PBEstruct: struct of PBE results
%       Fields: events, event_lengths, spike_order, frac_spike, ranks_vec,
%
% Optional:
%   - plotAll
%
% Detection method:
% 1) Detect periods of z-scored population activity rates above modelParam.PBE_zscore
%   standard deviations AND above a minimum rate of modelParam.PBE_min_Hz
% Optional) if modelParam.PBE_useMeanCrossing, move out the start/end of
%   the event to the mean-crossing points about the high-activity period
% 2) Combine candidate periods if they are within
%   modelParam.PBE_max_combine seconds of each other
% 4) Remove events that are shorter than modelParam.PBE_min_dur seconds
%


%% Parse inputs
inputObj = inputParser;
addRequired(inputObj, 'spikes',    @ismatrix )
addRequired(inputObj, 'modelParam',	@isstruct)
addParameter(inputObj, 'plotAll', false, @islogical);
parse(inputObj,spikes, modelParam, varargin{:});
p = inputObj.Results;

if size(spikes, 1)~=modelParam.n_E
    warning(['detect_PBE uses all cells in the spikes matrix. Did you mean', ...
        ' to pass in only E-cells?'])
end


%% Set up
tVec = 0:modelParam.dt:modelParam.t_max_preplay;    % Preplay time vector
popRateVec = smoothdata( mean(spikes, 1)/modelParam.dt, 'gaussian', modelParam.PBE_window); % Smoothed population firing rate firing rate
popRateMean = mean(popRateVec);
popRateSTD = std(popRateVec);
PBEthresh = max( popRateMean+(modelParam.PBE_zscore*popRateSTD), modelParam.PBE_min_Hz); % Threshold rate for PBE detection

if p.plotAll
    figure; plotSpikeRaster(spikes, 'TimePerBin', modelParam.dt, 'PlotType', 'scatter');
end


%% Detect PBE timepoints

% 1) Detect periods of high population activity
PBE_candidate1 = [popRateVec>PBEthresh];	% Binary vector of candidate PBE times
if p.plotAll
    [~, numRegions1] = bwlabel(PBE_candidate1);
    figure; plot(tVec, popRateVec); hold on; yline(popRateMean); yline(PBEthresh, 'r'); yyaxis right; plot(tVec, PBE_candidate1)
    title(['Step 1: ', num2str(numRegions1), ' events'])
    xlim([min(tVec), max(tVec)]); ylim([0, 1.05])
end

% Optional step: Set start/end at return to mean
if modelParam.PBE_useMeanCrossing
    PBE_candidateMeanCrossing = PBE_candidate1;
    onsetinds = find( diff(PBE_candidate1)==1)+1;
    while ~isempty(onsetinds)
        for i = numel(onsetinds):-1:1
            if popRateVec(onsetinds(i))>popRateMean
                onsetinds(i) = onsetinds(i)-1;
                PBE_candidateMeanCrossing(onsetinds(i)) = 1;
            else
                onsetinds(i) = [];
            end
        end
    end
    offsetinds = find( diff(PBE_candidate1)==-1);
    while ~isempty(offsetinds)
        for i = numel(offsetinds):-1:1
            if popRateVec(offsetinds(i))>popRateMean
                offsetinds(i) = offsetinds(i)+1;
                PBE_candidateMeanCrossing(offsetinds(i)) = 1;
            else
                offsetinds(i) = [];
            end
        end
    end
    if p.plotAll
        [~, numRegionsOptional] = bwlabel(PBE_candidateMeanCrossing);
        figure; plot(tVec, popRateVec); hold on; yline(popRateMean); yline(PBEthresh, 'r'); yyaxis right; plot(tVec, PBE_candidateMeanCrossing)
        title(['Optional, mean crossing thresh: ', num2str(numRegionsOptional), ' events'])
        xlim([min(tVec), max(tVec)]); ylim([0, 1.05])
    end
    PBE_candidate2 = PBE_candidateMeanCrossing;
else
    PBE_candidate2 = PBE_candidate1;
end

% 2) Combine candidate PBE that are very close in time
onsetInds = find(diff(PBE_candidate2)==1)+1;    % Time indices where a candidate PBE starts
offsetInds = find(diff(PBE_candidate2)==-1)+1;  % Time indices where a candidate PBE ends
if PBE_candidate2(1)==1; offsetInds = offsetInds(2:end); end    % If trial starts with PBE, ignore it
if PBE_candidate2(end)==1; onsetInds = onsetInds(1:end-1); end  % If trial ends with PBE, ignore it
for i = numel(offsetInds)-1:-1:1
    if onsetInds(i+1)-offsetInds(i)<modelParam.PBE_max_combine/modelParam.dt
        PBE_candidate2(offsetInds(i)-1:onsetInds(i+1)+1) = 1;
    end
end
if p.plotAll
    [~, numRegions2] = bwlabel(PBE_candidate2);
    figure; plot(tVec, popRateVec); hold on; yline(popRateMean); yline(PBEthresh, 'r'); yyaxis right; plot(tVec, PBE_candidate2)
    title(['Step 2: ', num2str(numRegions2), ' events'])
    xlim([min(tVec), max(tVec)]); ylim([0, 1.05])

    figure; histogram((offsetInds-onsetInds)*modelParam.dt);
    title('Event durations before removal of short events', 'FontWeight', 'Normal')
end

% 3) Remove candidate PBEs that don't last long enough
PBE_candidate3 = PBE_candidate2;
onsetInds = find(diff(PBE_candidate3)==1)+1; % indexes where a candidate PBE starts
for i = numel(onsetInds):-1:1
    if onsetInds(i)+modelParam.PBE_min_dur/modelParam.dt > numel(tVec)
        PBE_candidate3(onsetInds(i):end) = 0;
    elseif ~all(PBE_candidate3(onsetInds(i):onsetInds(i)+modelParam.PBE_min_dur/modelParam.dt))
        if i~=numel(onsetInds)
            PBE_candidate3(onsetInds(i):onsetInds(i+1)-1) = 0;
        else
            PBE_candidate3(onsetInds(i):onsetInds(i)+modelParam.PBE_min_dur/modelParam.dt) = 0;
        end
    end
end
if p.plotAll
    [~, numRegions3] = bwlabel(PBE_candidate3);
    figure; plot(tVec, popRateVec); hold on; yline(popRateMean); yline(PBEthresh, 'r'); yyaxis right; plot(tVec, PBE_candidate3)
    title(['Step 3: ', num2str(numRegions3), ' events'])
    xlim([min(tVec), max(tVec)]); ylim([0, 1.05])
end


%% Analyze and store PBE contents

onsetInds_final = find(diff(PBE_candidate3)==1)+1;	% Indexes where processed PBEs starts
offsetInds_final = find(diff(PBE_candidate3)==-1)+1; % Indexes where processed PBEs ends
if PBE_candidate3(1)==1; offsetInds_final = offsetInds_final(2:end); end    % If trial starts with PBE, ignore it
if PBE_candidate3(end)==1; onsetInds_final = onsetInds_final(1:end-1); end  % If trial ends with PBE, ignore it

if ~isempty(offsetInds_final) && offsetInds_final(end)>numel(popRateVec)
    warning(['Bug found. offsetInds end is ', num2str(offsetInds_final(end))])
    offsetInds_final(end) = numel(popRateVec);
end

events = [onsetInds_final', offsetInds_final'];
event_lengths = [offsetInds_final' - onsetInds_final']*modelParam.dt;

if size(events, 1)==0 % If no events, set results as empty
    PBEstruct.events = [];
    PBEstruct.event_lengths = [];
    PBEstruct.spike_order = [];
    PBEstruct.frac_spike = {};
    PBEstruct.ranks_vec = [];

else % If PBEs detected, calculate and store results
    PBEstruct.events = events;
    PBEstruct.event_lengths = event_lengths;

    %Find spike sequence for each event
    for ithEvent = 1:size(events, 1)
        event_spikes = spikes(:,events(ithEvent,1):events(ithEvent,2));
        % figure; plotSpikeRaster(event_spikes, 'TimePerBin', parameters.dt, 'PlotType', 'scatter');
        [e_spikes_x, ~] = find(event_spikes); % find neuron index of all spikes in order of spike-time
        spike_order = unique(e_spikes_x,'stable'); % take only first spike of each neuron
        ranks_vec = nan(1, size(spikes, 1));
        for k = 1:length(spike_order)
            n_ind = spike_order(k);
            ranks_vec(1,n_ind) = k;
        end
        PBEstruct.spike_order = spike_order;
        PBEstruct.frac_spike{ithEvent} = sum(~isnan(ranks_vec))/numel(ranks_vec);
        PBEstruct.ranks_vec(:,ithEvent) = ranks_vec;
    end
end


end