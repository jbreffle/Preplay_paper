function Config = setConfig()
% Config = utils.setConfig();
%
% Creates struct with all project-level parameters
% Loads data paths from config.json
%
% Paths for the project should be defined only in setConfig() if they are
% relative paths, and in config.json if they are machine-specific paths.
%
% Relatives paths are defined by the project root, which is defined by the
% location of the .git folder.


%% User directory

userDir = char(java.lang.System.getProperty('user.home'));


%% Local project paths

% git tracked local folders
projectPath = utils.getProjectPath();
functionsPath = fullfile(projectPath, 'functions');

% .gitignored local files, which might need to be created
gridSimJsonFilePath = fullfile(projectPath, 'gridSimNames.json');
if ~exist(gridSimJsonFilePath, 'file')
    copyfile('gridSimNames.json.example', 'gridSimNames.json')
end

% .gitignored local folders, which might need to be created
decodeDataPath = fullfile(projectPath, 'data', 'decodes');
if ~exist(decodeDataPath, 'dir')
    mkdir(decodeDataPath)
end


%% Read user configured raw and processed data paths from config.json

% Load config.json to struct
userConfigFilePath = fullfile(projectPath, 'config.json');
if exist(userConfigFilePath, 'file')
    userConfigFile = fileread(userConfigFilePath);
else
    error(['utils.setConfig() could not find the file config.json. ' ...
        'Did you create it following the example of ' ...
        'config.json.example?']);
end
userConfigStruct = jsondecode(userConfigFile);

% Increment through all paths in userConfigStruct and check if they exist
fields = fieldnames(userConfigStruct);
for i = 1:length(fields)
    % Is path if ends with "Folder"
    isFolder = endsWith(fields{i}, 'Path');
    if isFolder && ~exist(userConfigStruct.(fields{i}), 'dir')
        error(['utils.setConfig() could not find ' fields{i} '. ' ...
            'Did you set it in config.json?']);
    end
end


%% Set Config struct

Config = struct();
% local paths
Config.userDir = userDir;
Config.projectPath = projectPath;
Config.functionsPath = functionsPath;
Config.decodeDataPath = decodeDataPath;
% config.json
Config.simDataPath = fullfile(userConfigStruct.simDataPath);
Config.oldSimDataPath = fullfile(userConfigStruct.oldSimDataPath);
Config.exptDataPath = fullfile(userConfigStruct.exptDataPath);
Config.exptReplayDecodePath = fullfile(userConfigStruct.exptReplayDecodePath);
Config.exptPreplayDecodePath = fullfile(userConfigStruct.exptPreplayDecodePath);
Config.animalPrefixes = userConfigStruct.animalPrefixes';
Config.gridSimJsonFilePath = gridSimJsonFilePath;


%% Add project functions folder to path

addpath(fullfile(projectPath, 'src'))
addpath(fullfile(projectPath, 'src', 'jl-src'))


%% Confirm necessary add-ons

%TODO: define required add-ons
addons = matlab.addons.installedAddons;
%if ~ismember("SG", addons{:,"Identifier"})
%    warning('Signal Processing Toolbox is needed')
%end


end