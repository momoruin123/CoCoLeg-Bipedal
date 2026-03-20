function p_idx = getParameterIndex(paramName, config)
% GETPARAMETERINDEX Get parameter indices from model parameter vector
%
%   Inputs:
%     paramName - parameter name(s) as string or cell array of strings
%     config    - system configuration structure
%
%   Output:
%     p_idx     - indices of specified parameters in parameter vector
%
%  author: Iskandar Khemakhem, Oussama Barhoumi, IAMS, Uni Stuttgart, 2025

    % Get parameter names from model-specific function
    get_Parameters_names = str2func([config.model_name, '.getParameterNames']);
    parameterOrder = get_Parameters_names(config);

    % Initialize output
    if ischar(paramName) || isstring(paramName)
        paramName = {char(paramName)};  % convert to cell for uniformity
    end

    % Preallocate output
    p_idx = nan(size(paramName));

    for i = 1:numel(paramName)
        idx = find(strcmpi(parameterOrder, paramName{i}));
        if isempty(idx)
             warning('Parameter "%s" not found in standard parameter vector.', paramName{i});
             p_idx(i) = [];
        else
             p_idx(i) = idx;
        end
    end
end
