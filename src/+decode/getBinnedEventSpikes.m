function spkPerBin = getBinnedEventSpikes(celldata, riptimes, eventindex, event, modelParam)
% spkPerBin = decode.getBinnedEventSpikes();
%
%

%% Set up

sToMsConversion = 1000;
tBinSz = modelParam.tBinSz;
hpnum = numel(unique((celldata(:,3))));

%%
cellsiAll = celldata(find(celldata(:,2)==eventindex(event)),3);
[cellsi, ia] = unique(cellsiAll, 'first');

% Spikes per cell per time bin
%-----create the event matrix during SWRs (spkT{cells}.spiketimes) -----%
% Output: spkPerBin, nSpkPerTBin
[~,sortorder] = sort(ia);
event_cellSeq = cellsi(sortorder);
tmpind = find(celldata(:,2) == eventindex(event));
spiketimes = celldata(tmpind,1);
cellindex = celldata(tmpind,3);
for cell = event_cellSeq'
    validspikeidx = find(cellindex == cell);
    spkT{cell} = spiketimes(validspikeidx) * sToMsConversion;
end
startevent = riptimes(eventindex(event),1) * sToMsConversion;
endevent = riptimes(eventindex(event),2) * sToMsConversion;
timebins = startevent:tBinSz:endevent; % timebins are the binedges
nTBin = length(timebins)-1;
nCell = hpnum;
spkPerBin = zeros(1,nTBin, nCell); % keep the inactive cells as 0s.
for nn  = 1:hpnum
    cellInd = nn; %current cell
    if length(spkT) >= cellInd
        if ~isempty(spkT{cellInd})
            temp = histc(spkT{cellInd}, timebins); %[1 x nTBin x nCell]
            spkPerBin(1,:,cellInd) = temp(1:end-1);
        end
    end
end


end