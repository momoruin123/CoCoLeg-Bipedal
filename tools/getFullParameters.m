function [fullParams, paramNames] = getFullParameters(config)
% GETFULLPARAMETERS Generate complete parameter vector for specified model
%
%   Inputs:
%     config - system configuration structure
%
%   Outputs:
%     fullParams - vector of physical parameter values
%     paramNames - cell array of corresponding parameter names
%
%  author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

     if startsWith(config.model_name, 'h', 'IgnoreCase', true)
        [fullParams, paramNames] = getFullParameters_hopper(config);
     elseif startsWith(config.model_name, 'b', 'IgnoreCase', true)
        [fullParams, paramNames] = getFullParameters_biped(config);
     else
        error('Class not supported yet. Please adapt code to accomodate your class');
     end
end