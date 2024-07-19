function [ctxidx, hpidx] = getCellIdx(animalDir, animalprefix, sectionDay, epoch, isSimulatedData)
% [ctxidx, hpidx] = decode.getCellIdx(animalDir, animalprefix, sectionDay, epoch, isSimulatedData);
%
% Each row is formated as (tet, cell)

%% Set up

singleDay2019Animals = {'ER1','KL8','JS14','JS15','JS17','JS21'};


%% Calculate outputs

excludeList = decode.getExcludeList(animalprefix);

if isSimulatedData
    tetinfo = load(fullfile(sprintf('%s%stetinfo.mat',animalDir,animalprefix)), 'tetinfo'); % get num tets, get num cells
    tetinfo = tetinfo.tetinfo;
    numTets = numel(tetinfo{sectionDay}{epoch});
    numCells = tetinfo{sectionDay}{epoch}{1}.numcells;
    hpidx = [repmat(numTets, numCells, 1), [1:numCells]'];
    ctxidx = [];
elseif any(strcmp(animalprefix, singleDay2019Animals))
    % Function from Shin et al., 2019, for continuously tracked cell idx
    [ctxidx, hpidx] = matchidx_acrossep_singleday(animalDir, animalprefix, sectionDay, excludeList);
else
    % Find all cells for the given day and epoch
    cellinfo = loaddatastruct(animalDir,animalprefix,'cellinfo');
    nTets = numel(cellinfo{sectionDay}{epoch});
    hpidx = [];
    ctxidx = [];
    for ithTet = 1:nTets
        nCells = numel(cellinfo{sectionDay}{epoch}{ithTet});
        for ithCell = 1:nCells
            if ~isfield(cellinfo{sectionDay}{epoch}{ithTet}{ithCell}, 'tag2')
                continue
            else
                isCA1Pyr = isequal(cellinfo{sectionDay}{epoch}{ithTet}{ithCell}.tag2, 'CA1Pyr');
                if isCA1Pyr
                    hpidx = [hpidx; ithTet, ithCell];
                end
                isPFC = isequal(cellinfo{sectionDay}{epoch}{ithTet}{ithCell}.tag2, 'PFC');
                if isPFC
                    ctxidx = [ctxidx; ithTet, ithCell];
                end
            end
        end
    end
    % Remove rows that match exclude_list in either ctxidx or hpidx
    if ~isempty(excludeList)
        ctxidx = setdiff(ctxidx, excludeList,'rows');
        hpidx = setdiff(hpidx, excludeList,'rows');
    end
end
