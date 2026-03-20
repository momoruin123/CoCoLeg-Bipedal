function params = getModelParameters(config, casadiFlag, optFlag)
% GETMODELPARAMETERS Retrieve model parameters as symbolic variables or numerical values
%
%   Inputs:
%     config    - system configuration structure
%     casadiFlag- flag to return symbolic CasADi variables (true) or numeric values (false)
%     optFlag   - flag to exclude optimization parameters from the list
%
%   Output:
%     params    - vector of parameters (symbolic or numeric)
%  
%   author: Iskandar Khemakhem, Oussama Barhoumi, IAMS, Uni Stuttgart, 2025

    import casadi.* 

    get_Parameters_names = str2func([config.model_name, '.getParameterNames']);
    parameterList = get_Parameters_names(config); 

    %Eliminate optimization variables
    if optFlag
        %Get the list of the optimization parameters from config
        optParams = config.optParameterNames;

        if ~isempty(optParams)
            for i = 1:length(optParams)
                currentParam = optParams{i};
    
                % Find index of current parameter in parameterList
                idx = find(strcmp(parameterList, currentParam));
    
                % If found, remove from parameterList
                if ~isempty(idx)
                    parameterList(idx) = [];
                end
            end
        end
    end


    if casadiFlag
        % Initialize an empty cell array to store symbolic variables
        params = cell(length(parameterList), 1);

        % Iterate through each parameter name and create symbolic variables
        for i = 1:length(parameterList)
            paramName = parameterList{i};
            params{i} = MX.sym(paramName);  % Create symbolic variable
        end

        % Convert cell array to symbolic vector if needed
        params = [params{:}].';  % Transpose to get column vector
        
    else
        % Initialize an empty array to store symbolic variables
        params = zeros(length(parameterList), 1);

        %Get full parameters
        [fullParams, fullParamNames] = getFullParameters(config);

        % Check if config defines specific parameters and values
        hasOptParams = isfield(config, 'optParameterNames') && isfield(config, 'optParameterInit');

        % Iterate through each parameter in parameterList
        for i = 1:length(parameterList)
            currentParam = parameterList{i};
            if hasOptParams && any(strcmp(config.optParameterNames, currentParam))      
                % Find index in optParameterNames
                idx_opt = find(strcmp(config.optParameterNames, currentParam));
                params(i) = config.optParameterInit(idx_opt);
            else
                % Find the index of the current parameter in fullParamNames
                idx = find(strcmp(fullParamNames, currentParam));
            
                if isempty(idx)
                    error('Parameter "%s" not found in full parameter list', currentParam);
                end
            
                % Store the corresponding value
                params(i) = fullParams(idx);
            end
        end
    end
end

