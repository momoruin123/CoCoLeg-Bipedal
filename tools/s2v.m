function [vector, names] = s2v(struct, varargin)
% S2V Convert structure to vector with parameter names
%
%   Converts a structure of parameters into a numerical vector while
%   preserving parameter names. Handles both scalar and matrix-valued
%   parameters with automatic indexing.
%
%   Inputs:
%     struct   - parameter structure with numeric fields
%     names    - optional cell array of specific field names to extract
%
%   Outputs:
%     vector   - numerical vector of parameter values
%     names    - cell array of corresponding parameter names
%
%   Usage:
%     [v, n] = s2v(paramStruct)          % All parameters
%     [v, n] = s2v(paramStruct, {'m', 'g'}) % Specific parameters
%
% author: Korbinian Griesbauer, TUM.

    % Handle case where specific field names are provided
    if nargin > 1
        names = varargin{1};
        notFoundIndex = [];
        vector = zeros(length(names),1);
        for i = 1:length(names)
            try
                vector(i) = eval(['struct.',names{i},';']);
            catch ME
                if strcmp(ME.identifier,'MATLAB:nonExistentField')
                    notFoundIndex = [notFoundIndex, i];
                    vector(i) = 0;
                end
            end
        end
    else
        % Extract all fields from structure
        namesCell = fieldnames(struct);
        counter = 0;
        vector = [];
        for i = 1:length(namesCell)
            value = struct.(namesCell{i});
            [m, n] = size(value);
            % Handle scalar and matrix parameters
            for j = 1:m
                for k = 1:n
                    counter = counter+1;
                    vector(counter) = value(j,k);
                    if m*n == 1 % scalar parameter
                        names{counter} = namesCell{i};
                    else        % matrix parameter with indices
                        names{counter} = [namesCell{i},'(',num2str(j),',',num2str(k),')'];
                    end
                end
            end
        end
        vector = vector';
    end
end