function writeToGridSimNameJson(newSimFilename)
% writeToGridSimNameJson("grid_yyyy-mm-ddThh-MM")
%
% Adds a row to the JSON file gridSimName.json
% The row key is gridX where X is the ith row of the file


%%

% Set up
%if ~exist('+utils','dir'); addpath('..src'); end
Config = utils.setConfig;
gridSimJsonFilePath = Config.gridSimJsonFilePath;
rowIdentifierPrefix = 'grid';

% Load existing json
jsonText = fileread(gridSimJsonFilePath);
jsonData = jsondecode(jsonText);

% Add new row to json
jsonFieldNames = fieldnames(jsonData);
% Pattern: rowIdentifierPrefix<X>, where X is an integer, parse the jsonFieldNames to find the max X
maxGridInteger = 0;
for i = 1:length(jsonFieldNames)
    fieldName = jsonFieldNames{i};
    if startsWith(fieldName, rowIdentifierPrefix)
        % Extract the integer
        gridInteger = str2double(fieldName(length(rowIdentifierPrefix)+1:end));
        if gridInteger > maxGridInteger
            maxGridInteger = gridInteger;
        end
    end
end
nextFieldInteger = maxGridInteger+1;
newRowIdentifier = [rowIdentifierPrefix, num2str(nextFieldInteger)];
jsonData.(newRowIdentifier) = newSimFilename;

% Save modified json
jsonTextUpdated = jsonencode(jsonData, "PrettyPrint", true);
fid = fopen(gridSimJsonFilePath, 'w');
fprintf(fid, '%s', jsonTextUpdated);
fclose(fid);


end
