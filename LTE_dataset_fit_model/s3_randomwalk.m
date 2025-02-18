clc; clear; close all;

% Define dataset folders and class name

folders = {'./dataset/bus/', './dataset/car/'};
class = 'vehicle';

%folders = {'./dataset/static/'};
%class = 'static';

speed_limit = 20; % Speed threshold for filtering

% Initialize counters
total_transitions = 0; % Total number of valid CQI transitions
same_cqi_count = 0; % Number of times CQI stays the same

% Loop through each dataset folder
for f = 1:length(folders)
    dataset_dir = folders{f};

    % Get list of CSV files
    filePattern = fullfile(dataset_dir, '*.csv'); 
    csvFiles = dir(filePattern);
    
    % Process each CSV file
    for k = 1:length(csvFiles)
        filePath = fullfile(csvFiles(k).folder, csvFiles(k).name);
        
        try
            % Read CSV file
            data = readtable(filePath, 'TextType', 'string');
            
            % Ensure required columns exist
            if all(ismember({'NetworkMode', 'CQI', 'Speed'}, data.Properties.VariableNames))
                % Convert NetworkMode to string
                if ~isa(data.NetworkMode, 'string')
                    data.NetworkMode = string(data.NetworkMode);
                end
                
                % Convert CQI to numeric, replacing "-" with NaN
                if isa(data.CQI, 'string')
                    data.CQI(data.CQI == "-") = missing;
                    data.CQI = str2double(data.CQI);
                end
                
                % Convert Speed to numeric
                if isa(data.Speed, 'string')
                    data.Speed = str2double(data.Speed);
                end
                
                % Filter LTE-only rows and Speed >= speed_limit
                validRows = (data.NetworkMode == "LTE") & (data.Speed >= speed_limit) & ~isnan(data.CQI);
                filteredData = data(validRows, :);
                
                % Extract CQI sequence
                cqi_sequence = filteredData.CQI;
                
                % Compute transitions
                for i = 1:length(cqi_sequence) - 1
                    if cqi_sequence(i+1) == cqi_sequence(i)
                        same_cqi_count = same_cqi_count + 1;
                    end
                    total_transitions = total_transitions + 1;
                end
            end
        catch ME
            fprintf('Error reading file %s: %s\n', filePath, ME.message);
        end
    end
end

% Compute probability
p_same_cqi = same_cqi_count / total_transitions;

% Display probability
fprintf('\nProbability of Staying on the Same CQI Value:\n');
fprintf('---------------------------------------------\n');
fprintf('P(CQI remains the same) = %.4f\n', p_same_cqi);

% Save results accordingly
if strcmp(class, 'static')
    p_same_cqi_static = p_same_cqi;
    save('cqi_stay_probability_static.mat', 'p_same_cqi_static');
    fprintf('\nProbability saved to "cqi_stay_probability_static.mat".\n');
else
    p_same_cqi_vehicle = p_same_cqi;
    save('cqi_stay_probability_vehicle.mat', 'p_same_cqi_vehicle');
    fprintf('\nProbability saved to "cqi_stay_probability_vehicle.mat".\n');
end
