function projectPath = getProjectPath()
% projectPath = utils.getProjectPath();
%
% Returns the path to the project folder
% Iteratively go up the directory tree until the .git folder is found,
% which defines the project root


%%  Get current directory
currDir = mfilename('fullpath');

%% Iterate up the directory tree until .git folder is found

while ~exist(fullfile(currDir, '.git'), 'dir')
    % Go up one directory
    currDir = fileparts(currDir);
    % If we've reached the computer's root directory, throw error
    if strcmp(currDir, filesep)
        error('Could not find project root directory with .git folder');
    end
end


%% currDir is the parent folder with the .git folder

projectPath = currDir;


end
