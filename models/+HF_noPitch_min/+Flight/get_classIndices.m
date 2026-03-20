function [idx_x, idx_u] = get_classIndices()
%GET_CLASSINDICES Returns the indices of the optimal states and inputs in the full states

% State indices
idx_x = [1; 2; 4; 5; 6; 7; 9; 10];

% Input indices
idx_u = [1; 2];

end
