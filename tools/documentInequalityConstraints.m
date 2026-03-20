function documentInequalityConstraints(config, outputFileName)
% EXPORTINEQUALITYCONSTRAINTS Write inequality constraint functions to a text file.
%
%   exportInequalityConstraints(config, outputFileName)
%
%   Inputs:
%       config          - struct with fields:
%                           model_name (string)
%                           phaseSequence (cell array of strings)
%                           useInequalityConstraints (logical)
%       outputFileName  - name of the output text file (string)
%
%  author: Iskandar Khemakhem, IAMS, Uni Stuttgart, 2025


    modelName     = config.model_name;
    phaseSequence = config.phaseSequence;
    useIneq       = config.useInequalityConstraints;

    fid = fopen(outputFileName, 'w');
    if fid == -1
        error('Could not open file %s for writing.', outputFileName);
    end
    
    if ~useIneq
        fprintf(fid, 'Inequality constraints are disabled (useIneq = false). No file written.\n');
        fclose(fid);
        return;
    end


    try
        for i = 1:numel(phaseSequence)
            phase = phaseSequence{i};
            modelPhaseName = [modelName, '.', phase];
            funcFullName = [modelPhaseName, '.inequaliyConstraints'];
    
            % Get the function file path using `which`
            funcPath = which(funcFullName);
    
            fprintf(fid, '########## Phase: %s ##########\n', phase);
    
            if ~isempty(funcPath)
                % Read and write function content
                fileContent = fileread(funcPath);
                fprintf(fid, '%s\n\n', fileContent);
            else
                fprintf(fid, 'No file found for %s\n\n', funcFullName);
            end
        end
    catch ME
        fclose(fid);
        rethrow(ME);
    end


    fclose(fid);
    fprintf('Inequality constraints written to %s\n', outputFileName);
end
