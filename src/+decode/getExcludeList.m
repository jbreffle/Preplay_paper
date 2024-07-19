function excludeList = getExcludeList(animalprefix)
% excludeList = decode.getExcludeList(animalprefix);
%
% List of interneurons (from Wenbo's code for Shin et al., 2019)

if strcmp(animalprefix,'ER1')
    excludeList = [1,1;1,2];
elseif strcmp(animalprefix,'KL8')
    excludeList = [];
elseif strcmp(animalprefix,'JS14')
    excludeList = [6,1;8,3;23,1];
elseif strcmp(animalprefix,'JS15')
    excludeList = [5,4;8,6];
elseif strcmp(animalprefix,'JS17')
    excludeList = [6,2;6,4;6,6;7,5;10,1;11,2;23,2;23,3;23,4];
elseif strcmp(animalprefix,'JS21')
    excludeList = [6,5;21,2;25,2;25,3];
else
    excludeList = [];
end
