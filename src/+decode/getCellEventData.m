function [celldata, eventindex] = getCellEventData(day, epoch, modelParam, riptimes, cellidxm, spikes)
% [celldata, eventindex] = decode.getCellEventData(day, epoch, modelParam, riptimes, cellidxm, spikes)
%
%

%% Set up

cellcountthresh = modelParam.cellcountthresh;
nCells = size(cellidxm, 1);


%% Return empty if no valid ripples

if isempty(riptimes)
    celldata = [];
    eventindex = [];
    return
end


%% Accumulate celldata and calculate active cells per event

celldata = [];
spikecounts = [];
for cellcount = 1:nCells
    index = [day,epoch,cellidxm(cellcount,:)] ;
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}) && ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
    else
        spiketimes = [];
    end
    spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
    if ~isempty(spiketimes)
        validspikes = find(spikebins);
        spiketimes = spiketimes(validspikes);
        spikebins = spikebins(validspikes);
        tmpcelldata = [spiketimes spikebins];
    end
    if ~isempty(spiketimes)
        tmpcelldata(:,3) = cellcount;
    else
        tmpcelldata = [0 0 cellcount];
    end
    celldata = [celldata; tmpcelldata];
    spikecount = zeros(1,size(riptimes,1));
    for i = 1:length(spikebins)
        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
    end
    spikecounts = [spikecounts; spikecount];
end
cellcounts = sum((spikecounts > 0));
eventindex = find(cellcounts >= cellcountthresh); % exclude events < cellcountthresh active


end