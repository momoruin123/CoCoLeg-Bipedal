function [succStruct, z_final, status, cost] = continuation_strategy(grid_array, cost, status, indices, costRun, z_finalRun, convergedRun, z_final, n_batch, config)
%CONTINUATION_STRATEGY Router function for continuation strategies
%
% This function serves as a dispatcher that selects and executes the
% appropriate continuation strategy based on configuration settings.
%
%   Inputs:
%     grid_array   - M×N matrix of grid points (M points, N parameters)
%     cost         - M×1 vector of cost values for each grid point
%     status       - M×1 vector of status codes for each grid point
%     indices      - indices of current batch being processed
%     costRun      - cost values from current optimization batch
%     z_finalRun   - decision variables from current optimization batch
%     convergedRun - convergence flags from current optimization batch
%     z_final      - M×n_z matrix of final decision variables
%     n_batch      - number of points in current batch
%     config       - configuration structure with continuation parameters
%
%   Outputs:
%     succStruct   - structure containing successor information:
%                    status: status codes for successors (-1=similar exists, 0=new)
%                    gridPts: grid coordinates of successors
%                    cost: cost values from parent solutions
%                    z_init: initial guesses for successors (parent solutions)
%                    z_final: placeholder for successor solutions
%     z_final      - updated decision variable matrix
%     status       - updated status vector
%     cost         - updated cost vector
%   Status Codes:
%     0: Untried   1: Optimal   2: Suboptimal   3: Failed
%     -1: Successor with similar solution exists
%
% STRATEGY SELECTION:
% Currently only 'poorman' strategy is implemented, which processes
% points sequentially using neighboring solutions as initial guesses.
%
%  author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

% Select and execute the appropriate continuation strategy
if strcmpi(config.cont.type, 'poorman')
    [succStruct, z_final, status, cost] = continuation_strategy_poorman(grid_array, cost, status, indices, costRun, z_finalRun, convergedRun, z_final, n_batch, config);
else
    error('Continuation strategy not supported yet. Please adapt code to accomodate your strategy');
end

end