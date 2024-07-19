function [rm, pm, tm, cellidxm] = createSpatialMatrices(animalprefix, hpidx, loadDir, day, eprun, modelParam)
% [rm, pm, tm, cellidxm] = decode.createSpatialMatrices(animalprefix, day, epoch);
%
%-----create the ratemaps [nPosBin x nHPCells]-----%
% rm: ratemap matrix
% pm: position matrix
% tm: track matrix
% cellidxm: cell index matrix

minPeakRate = modelParam.minPeakRate; % Hz, minimum peak rate to include cell as Place Cell
wellcutoff = modelParam.wellcutoff;         % cm, remove reward-well regions (15cm around start and end); or 0cm without exclusion

rm = []; % ratemap matrix
pm = []; % position matrix
tm = []; % track matrix
cellidxm = [];
%load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
linfields = load(sprintf('%s%slinfields0%d.mat',loadDir,animalprefix,day), 'linfields'); % get num tets, get num cells
linfields = linfields.linfields;

% nTraj = 4;
nTraj = numel(linfields{day}{eprun}{hpidx(1,1)}{hpidx(1,2)});

hpnum = length(hpidx(:,1));
for i = 1:hpnum
    cind = hpidx(i,:);
    if (length(linfields{day}{eprun})>= cind(1))
        if  (length(linfields{day}{eprun}{cind(1)})>= cind(2))
            linfield1 = linfields{day}{eprun}{cind(1)}{cind(2)};
        else
            linfield1 =[];
        end
    else
        linfield1=[];
    end

    if ~isempty(linfield1)
        linfield_hp = [];
        lintrack_hp = [];
        pos_hp = [];
        for track = 1:nTraj
            temp1 = linfield1{track};
            pos1 = temp1(:,1);
            lintrack1 = ones(size(pos1))*track;
            occnormrate1 = temp1(:,5);
            linfield_hp = [linfield_hp;occnormrate1];
            pos_hp = [pos_hp;pos1];
            lintrack_hp = [lintrack_hp;lintrack1];
        end
        if (max(linfield_hp) >= minPeakRate) % peak firing rate max larger than minPeakRate
            a = find(isnan(linfield_hp));
            %pad nan
            if ~isempty(a)
                [lo,hi]= findcontiguous(a);  %find contiguous NaNs
                for ii = 1:length(lo)
                    if lo(ii) > 1 & hi(ii) < length(linfield_hp)
                        fill = linspace(linfield_hp(lo(ii)-1), ...
                            linfield_hp(hi(ii)+1), hi(ii)-lo(ii)+1);
                        linfield_hp(lo(ii):hi(ii)) = fill;
                    end
                end
            end
            rm = [rm;linfield_hp'];
            pm = [pm;pos_hp'];
            tm = [tm;lintrack_hp'];
            cellidxm = [cellidxm; cind];
        end
    end
end
rm = rm'; %[nPosBin x nHPCells]
pm = pm';
tm = tm';

% remove reward-well regions (15cm around start and end)
for i = 1:nTraj
    pm_traj = pm(find(tm == i));
    maxpos = max(max(pm_traj));
    rm(find(tm == i & pm <= wellcutoff)) = 0;
    rm(find(tm == i & pm >= maxpos-wellcutoff)) = 0;
end

rm = rm+ (eps.^8); %Add a small number so there are no zeros

end