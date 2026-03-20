function varargout = save_data(basename, saveMode, varargin)
% SAVE_DATA Save data to .mat files with flexible storage modes
%
%   Inputs:
%     basename - base filename for saving
%     saveMode - 'overwrite', 'version', or 'append'
%     varargin - name-value pairs of variables to save
%
%   Outputs:
%     filename - full path to saved file (optional)
%
%   Features:
%     - Automatic organization in 'data' subfolder
%     - Version numbering for 'version' mode
%     - Conflict resolution with indexed suffixes for 'append' mode
%     - Metadata extraction from config structure
%
%   author: Iskandar Khemakhem, Oussama Barhoumi, IAMS, Uni Stuttgart, 2025

% Get path to data subfolder in parent directory
thisFilePath = mfilename('fullpath');
[thisFolder, ~, ~] = fileparts(thisFilePath);
parentFolder = fileparts(thisFolder);
dataFolder = fullfile(parentFolder, 'data');

try
    % Create data folder if it doesn't exist
    if ~isfolder(dataFolder)
        mkdir(dataFolder);
    end

    % Validate name-value pair formatting
    if mod(length(varargin), 2) ~= 0
        error('Arguments must be in name-value pairs.');
    end

    % Initialize metadata fields
    costName = '';
    costNorm = '';
    modelName = '';
    config = [];

    % Extract metadata from config structure if provided
    for i = 1:2:length(varargin)
        name = varargin{i};
        value = varargin{i+1};
        if strcmp(name, 'config')
            config = value;
        end
    end

    % Populate metadata from config structure
    if ~isempty(config)
        if isempty(costName) && isfield(config, 'costName')
            costName = config.costName;
        end
        if isempty(costNorm) && isfield(config, 'costNormalization')
            costNorm = config.costNormalization;
        end
        if isempty(modelName) && isfield(config, 'model_name')
            modelName = config.model_name;
        end
    end

    % Construct filename prefix from metadata
    parts = {basename};
    if ~isempty(costName), parts{end+1} = costName; end
    if ~isempty(costNorm), parts{end+1} = costNorm; end
    if ~isempty(modelName), parts{end+1} = modelName; end
    basePrefix = strjoin(parts, '_');
    
    % Find existing files matching the base pattern
    filePattern = fullfile(dataFolder, [basePrefix '*.mat']);
    existingFiles = dir(filePattern);
    
    % Determine filename based on save mode
    switch lower(saveMode(1))  % Check first character only
        case 'o'  % overwrite
            filename = fullfile(dataFolder, [basePrefix '.mat']);
            
        case 'v'  % versioning
            % Find highest existing version number
            maxVersion = 0;
            for i = 1:length(existingFiles)
                [~, name] = fileparts(existingFiles(i).name);
                versionPart = regexp(name, ['^' basePrefix '_(\d+)$'], 'tokens');
                if ~isempty(versionPart)
                    versionNum = str2double(versionPart{1}{1});
                    if versionNum > maxVersion
                        maxVersion = versionNum;
                    end
                end
            end
            filename = fullfile(dataFolder, sprintf('%s_%d.mat', basePrefix, maxVersion + 1));
            
        case 'a'  % append
            filename = fullfile(dataFolder, [basePrefix '.mat']);
            if isfile(filename)
                % Load existing variables and handle naming conflicts
                existingVars = load(filename);
                
                for i = 1:2:length(varargin)
                    varName = varargin{i};
                    varValue = varargin{i+1};
                    
                    if isfield(existingVars, varName)
                        % Find existing indexed versions of this variable
                        allVars = fieldnames(existingVars);
                        pattern = ['^' varName '__i(\d+)$'];
                        indices = [];
                        
                        % Check for base variable and indexed versions
                        for j = 1:length(allVars)
                            tokens = regexp(allVars{j}, pattern, 'tokens');
                            if ~isempty(tokens)
                                indices(end+1) = str2double(tokens{1}{1});
                            end
                        end
                        
                        % Determine next available index
                        if isempty(indices)
                            newIndex = 1;
                        else
                            newIndex = max(indices) + 1;
                        end
                        
                        % Create indexed variable name
                        newName = sprintf('%s__i%d', varName, newIndex);
                        existingVars.(newName) = varValue;
                    else
                        % No conflict - use original name
                        existingVars.(varName) = varValue;
                    end
                end
                
                % Save merged data and return
                save(filename, '-struct', 'existingVars', '-v7.3');
                
                if nargout >= 1
                    varargout{1} = filename;
                end
                return;
            end
            % Fall through to normal save if no existing file
            
        otherwise
            error('Invalid save mode. Use ''overwrite'', ''version'', or ''append''.');
    end

    % Standard save for non-append cases
    varNames = varargin(1:2:end);
    varValues = varargin(2:2:end);
    for i = 1:length(varNames)
        eval([varNames{i} ' = varValues{i};']);
    end
    
    save(filename, varNames{:}, '-v7.3');

    if nargout >= 1
        varargout{1} = filename;
    end
    
catch ME
    errordlg(['Could not save: ', ME.message], 'Save Error');
    if nargout >= 1
        varargout{1} = '';
    end
end
end