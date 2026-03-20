function paramNames = getParameterNames(config)
%GETPARAMETERNAMES_HOPPER Returns the parameter names for the Hopper model

paramNamesAll = {
    'g';
    'm';
    'l';
    'm_2';
    'm_3';
    'theta_1';
    'theta_2';
    'theta_3';
    'l_1';
    'l_2';
    'd_1';
    'd_2';
    'd_3';
    'alpha0';
    'beta0';
    'xi_h';
    'xi_k';
    'k_h';
    'k_k';
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
