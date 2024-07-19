function filename = getDecodeFilename(fileIndex)
% filename = getDecodeFile(fileIndex);
%


%% Load json file

Config = utils.setConfig;
gridSimJsonFilePath = Config.gridSimJsonFilePath;
jsonText = fileread(gridSimJsonFilePath);
jsonData = jsondecode(jsonText);


%% Get indexed row from jsonData is isfield(), else load from decodeFile

if isfield(jsonData, fileIndex)
    filename = jsonData.(fileIndex);
elseif isKey(decodeFile, fileIndex)
    filename = decodeFile(fileIndex);
    return
else
    error('File index not found.');
end


end