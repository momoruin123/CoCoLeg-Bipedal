function paramNames = getParameterNames(config)
%GETPARAMETERNAMES_HOPPER Returns the parameter names for the Hopper model
%    author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

paramNamesAll = {
    'm';
    'g';
    'l_0';
    'm_l';
    'm_f';
    'd_l';
    'd_f';
    'r_f';
    'xi_l';
    'k_l';
};

% Reorder if config.optParameterNames is given
if isfield(config, 'optParameterNames') && ~isempty(config.optParameterNames)
    optNames = config.optParameterNames;
    % Keep only those optimization names that exist
    optNames = optNames(ismember(optNames, paramNamesAll));
    % Remove optimization names from main list
    paramNames = setdiff(paramNamesAll, optNames, 'stable');
    % Append optimization names at the end
    paramNames = [paramNames; optNames(:)];
else
    paramNames = paramNamesAll;
end

end
