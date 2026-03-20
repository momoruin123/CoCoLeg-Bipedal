function extract_data(matFilePath)
% EXTRACT_DATA Load and combine fragmented variables from .mat file
%
%   Processes .mat files containing variables that have been split across
%   multiple fragments with '__i*' suffixes (e.g., data__i1, data__i2).
%   Combines fragmented variables along the first dimension and loads all
%   variables into the base workspace for analysis.
%
%   Inputs:
%     matFilePath - path to the .mat file (extension optional)
%
%   Behavior:
%     - Automatically appends .mat extension if not provided
%     - Identifies variables with '__i*' suffixes (e.g., data__i1, data__i2)
%     - Combines fragmented variables along first dimension using cat(1,...)
%     - Assigns combined result to base variable name (without suffix)
%     - Clears individual fragments from base workspace
%     - Loads all non-fragmented variables directly to base workspace
%
%   Example:
%     Given file 'results.mat' containing:
%        data__i1 = [1 2 3], data__i2 = [4 5 6], other_var = 42
%     After extract_data('results.mat'):
%        data = [1 2 3; 4 5 6], other_var = 42 in base workspace
%
%
% author: skandar Khemakhem, Oussema Barhoumi, IAMS, Uni Stuttgart, 2025

    % Ensure file has .mat extension
    if ~endsWith(matFilePath, '.mat')
        matFilePath = [matFilePath, '.mat'];
    end

    % Load all variables from the .mat file
    loadedData = load(matFilePath);
    vars = fieldnames(loadedData);

    % Find base variables that have indexed versions
    baseVars = {};
    for i = 1:length(vars)
        tokens = regexp(vars{i}, '^(.*?)__i\d+$', 'tokens');
        if ~isempty(tokens)
            baseVar = tokens{1}{1};
            if ~ismember(baseVar, baseVars)
                baseVars{end+1} = baseVar;
            end
        end
    end

    % Process each base variable with indexed fragments
    for i = 1:length(baseVars)
        baseVar = baseVars{i};

        % Find all versions of this variable (base + indexed)
        pattern = ['^' baseVar '(__i\d+)?$'];
        matchingVars = vars(~cellfun(@isempty, regexp(vars, pattern)));

        % Collect all arrays from fragments
        arrays = {};
        for j = 1:length(matchingVars)
            arrays{end+1} = loadedData.(matchingVars{j});
        end

        % Concatenate fragments along the first dimension
        combined = cat(1, arrays{:});

        % Assign combined result to base workspace
        assignin('base', baseVar, combined);

        % Clean up individual fragments from base workspace
        for j = 1:length(matchingVars)
            if ~strcmp(matchingVars{j}, baseVar)
                evalin('base', ['clear ', matchingVars{j}]);
            end
        end
    end

    % Assign non-indexed variables that weren't part of combination
    assignedVars = [baseVars, strcat(baseVars, '__i1')];  % Track assigned variables
    for i = 1:length(vars)
        varName = vars{i};
        % Skip indexed parts and already assigned base variables
        isIndexed = ~isempty(regexp(varName, '__i\d+$', 'once'));
        isBaseAssigned = ismember(varName, baseVars);
        if ~isIndexed && ~isBaseAssigned
            assignin('base', varName, loadedData.(varName));
        end
    end
end
