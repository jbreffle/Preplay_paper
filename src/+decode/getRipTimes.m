function riptimes = getRipTimes(animalprefix, loadDir, day, epoch, modelParam)
% riptimes = decode.getRipTimes(animalprefix, loadDir, day, epoch, modelParam);
%
%

%% Set up

minEventDur = modelParam.minEventDur;
sToMsConversion = 1000;


%% get ripple time

% load(sprintf('%s%srippletime0%d.mat',dir,animalprefix,day));

if isfield(modelParam, 'useMESrippletime') && modelParam.useMESrippletime
    ripple = load(sprintf('%s%srippletimes_MES.mat',loadDir,animalprefix), 'rippletimes_MES'); % get num tets, get num cells
    ripple = ripple.rippletimes_MES;
    disp(['Note: using ', sprintf('%srippletimes_MES.mat',animalprefix), ' file'])
else
    ripple = load(sprintf('%s%srippletime0%d.mat',loadDir,animalprefix,day), 'ripple'); % get num tets, get num cells
    ripple = ripple.ripple;
end
rip = ripple{day}{epoch};

if numel(rip.starttime)>numel(rip.endtime)
    riptimes(:,1) = rip.starttime(1:numel(rip.endtime));
    riptimes(:,2) = rip.endtime;
else
    riptimes(:,1) = rip.starttime;
    riptimes(:,2) = rip.endtime;
end

dur = (riptimes(:,2) - riptimes(:,1)) * sToMsConversion; % event duration in ms
keepidx = find(dur >= minEventDur); % exclude events < minEventDur ms
riptimes = riptimes(keepidx,:);


end