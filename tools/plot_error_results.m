function plot_error_results(config, solution, error, dt)
% PLOT_ERROR_RESULTS Visualize dynamics constraint violations and integration errors
%
%   Inputs:
%     config   - system configuration structure
%     solution - optimization solution variables
%     error    - error structure from analyze_error function
%     dt       - time step for trajectory interpolation
%
% author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025

    % Extract phase information from configuration
    nPhases = numel(config.phaseSequence);
    N_array = config.N;

    % Interpolate solution to fine trajectory with state derivatives
    traj = interpolate_Z2traj(config, solution, [], 1, dt, 1);

    % Identify all states that appear in any phase for comprehensive plotting
    max_n_x = numel(config.defaultStateNames);
    state_active_flags = false(1, max_n_x);
    
    % Mark states that are active in any phase
    for i = 1:nPhases
        phase_name = config.phaseSequence{i};
        getIndices = str2func([config.model_name, '.', phase_name, '.get_classIndices']);
        [idx_x, ~] = getIndices();
        state_active_flags(idx_x) = true;
    end

    % Get indices of all relevant states across all phases
    all_relevant_states = find(state_active_flags);
    n_total_states = numel(all_relevant_states);

    %% FIGURE 1: Instantaneous Dynamics Error vs Time
    figure('Name', 'Dynamics Error per State'); 
    for s_idx = 1:n_total_states
        j = all_relevant_states(s_idx);  % Global state index
        subplot(2, n_total_states/2, s_idx); 
        hold on; grid on;
        title(config.defaultStateNames{j});
        ylabel('Dynamic error'); 
        xlabel('Time');
    
        t_end = 0;  % Cumulative time offset for multi-phase plotting
        for i = 1:nPhases
            phase_name = config.phaseSequence{i};
            getIndices = str2func([config.model_name, '.', phase_name, '.get_classIndices']);
            [idx_x, ~] = getIndices();
    
            t_vec = traj(i).t;
            N_i = N_array(i);
            h = traj(i).t(end)/N_i;
            colloc_times = (0:N_i)*h + t_end;  % Collocation point times
    
            % Plot error if state exists in this phase
            if ismember(j, idx_x)
                e_vec = error(i).dyn_e(:, find(j == idx_x));
                plot(t_vec + t_end, e_vec, 'DisplayName', ['[', num2str(i), '] ', phase_name]);
                scatter(colloc_times, zeros(size(colloc_times)), 5, 'ko', 'HandleVisibility', 'off');
            else
                % Plot invisible line to maintain time axis consistency
                plot(t_vec + t_end, zeros(size(t_vec)), 'HandleVisibility','off', ...
                     'DisplayName', ['[', num2str(i), '] ', phase_name]);
            end
            t_end = t_end + t_vec(end);  % Update cumulative time
        end
        legend('show');
        hold off;
    end

    %% FIGURE 2: Romberg Integrated Error per Segment
    figure('Name', 'Romberg Integrated Error per State'); 
    for s_idx = 1:n_total_states
        j = all_relevant_states(s_idx);  % Global state index
        subplot(2, n_total_states/2, s_idx); 
        hold on; grid on;
        title(config.defaultStateNames{j});
        ylabel('Integrated Error'); 
        xlabel('Segment #');
    
        segment_counter = 0;  % Global segment counter across phases
    
        for i = 1:nPhases
            phase_name = config.phaseSequence{i};
            getIndices = str2func([config.model_name, '.', phase_name, '.get_classIndices']);
            [idx_x, ~] = getIndices();
    
            N_i = N_array(i);
            seg_indices = (1:N_i) + segment_counter;  % Global segment indices
    
            % Plot integrated error if state exists in this phase
            if ismember(j, idx_x)
                bar_vals = error(i).state_e(:, find(j == idx_x));
                scatter(seg_indices, bar_vals, 30, 'filled', ...
                        'DisplayName', ['[', num2str(i), '] ', phase_name]);
            else
                % Plot invisible points to maintain segment numbering
                scatter(seg_indices, zeros(size(seg_indices)), 1, 'filled', ...
                        'HandleVisibility', 'off');
            end
    
            % Mark phase boundaries with vertical lines
            xline(segment_counter + N_i + 0.5, 'r--', 'HandleVisibility', 'off');
            segment_counter = segment_counter + N_i;  % Update segment counter
        end
        xlim([0.5, segment_counter + 0.5]);  % Set x-axis limits
        legend('show');
        hold off;
    end
end