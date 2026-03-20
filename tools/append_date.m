function output_str = append_date(input_str)
% APPEND_DATE_TO_STRING Prepends the current date (day_month as numbers) to a string.
%   OUTPUT_STR = APPEND_DATE_TO_STRING(INPUT_STR) adds the current date in 
%   'DD_MM' format (e.g., '02_06') to the beginning of INPUT_STR.
%
%   Example:
%       append_date_to_string('logfile') → '02_06_logfile' (if run on June 2nd)
%
%  author: Oussema Barhoumi, IAMS, Uni Stuttgart, 2025

    % Get current date as day and month (numeric format)
    current_date = datetime('now', 'Format', 'dd_MM');
    date_str = char(current_date); % Convert to string (e.g., '02_06')
    
    % Append date to the input string with an underscore
    output_str = [date_str '_' input_str];
    output_str = output_str(:)';
end