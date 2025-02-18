clc; clear; close all;

% Define dataset directory
dataset_dir = './dataset';
output_dir = './nb_fits'; % New directory for Negative Binomial fits
if ~exist(output_dir, 'dir')
    mkdir(output_dir); % Create folder if it doesn't exist
end

% Get list of all CSV files inside all subdirectories
filePattern = fullfile(dataset_dir, '**', '*.csv'); % Search recursively
csvFiles = dir(filePattern);

% Initialize storage for extracted data
all_CQI = [];
all_DL_bitrate = [];

% Loop through each CSV file
for k = 1:length(csvFiles)
    filePath = fullfile(csvFiles(k).folder, csvFiles(k).name);
    
    % Read CSV file
    try
        data = readtable(filePath, 'TextType', 'string');

        % Ensure required columns exist
        if all(ismember({'NetworkMode', 'CQI', 'DL_bitrate'}, data.Properties.VariableNames))
            % Convert "NetworkMode" to string if needed
            if ~isa(data.NetworkMode, 'string')
                data.NetworkMode = string(data.NetworkMode);
            end
            
            % Convert "CQI" to numeric, replacing "-" or non-numeric values with NaN
            if isa(data.CQI, 'string')
                data.CQI(data.CQI == "-") = missing; % Convert "-" to missing
                data.CQI = str2double(data.CQI); % Convert to numeric
            end
            
            % Convert "DL_bitrate" to numeric
            if isa(data.DL_bitrate, 'string')
                data.DL_bitrate = str2double(data.DL_bitrate);
            end
            
            % Filter: Keep only LTE rows and valid CQI values
            validRows = (data.Operatorname == "A") & (data.NetworkMode == "LTE") & ~isnan(data.CQI);
            filteredData = data(validRows, :);
            
            % Store extracted values
            all_CQI = [all_CQI; filteredData.CQI];
            all_DL_bitrate = [all_DL_bitrate; filteredData.DL_bitrate];
        end
    catch ME
        fprintf('Error reading file %s: %s\n', filePath, ME.message);
    end
end

% Remove NaN values and convert bitrate to discrete counts
validIndices = ~isnan(all_CQI) & ~isnan(all_DL_bitrate);
all_CQI = all_CQI(validIndices);
all_DL_bitrate = round(all_DL_bitrate(validIndices) / 1e3); % Convert to kbps and round to discrete values

% Fit Negative Binomial distribution for each CQI level (0 to 15)
nbParams = cell(16, 2); % Store CQI and its NB parameters
fprintf('\nNegative Binomial Distribution Parameters (r, p) per CQI level:\n');
fprintf('-------------------------------------------------------------\n');

for cqi = 0:15
    % Extract DL_bitrate for this CQI level
    bitrate_cqi = all_DL_bitrate(all_CQI == cqi);
    
    if length(bitrate_cqi) > 10 % Ensure sufficient data points
        % Estimate Negative Binomial parameters using method of moments
        mean_bitrate = mean(bitrate_cqi);
        var_bitrate = var(bitrate_cqi);
        
        if var_bitrate > mean_bitrate % Validity condition for NB distribution
            r = mean_bitrate^2 / (var_bitrate - mean_bitrate); % Shape parameter
            p = mean_bitrate / var_bitrate; % Probability parameter

            % Store results
            nbParams{cqi+1, 1} = cqi; % CQI
            nbParams{cqi+1, 2} = [r, p]; % Negative Binomial parameters
            
            % Display result
            fprintf('CQI %2d: r = %.3f, p = %.3f\n', cqi, r, p);

            % Create plot
            fig = figure();
            ax = axes(fig);

            histogram(ax, bitrate_cqi, 'Normalization', 'pdf', 'FaceAlpha', 0.5);
            hold on;

            % Generate Gamma PDF
            x = [1:100];
            nb_pmf = nbinpdf(x, r, p);

            % Plot Gamma PDF
            plot(ax, x, nb_pmf, 'r-', 'LineWidth', 2);

            font = 'Times New Roman';
            set(fig,'defaultAxesFontName',font);
            set(fig,'DefaultTextFontName', font, 'DefaultAxesFontName', font);
            set(ax,'FontName', font, 'FontName', font);
            set(fig,'defaultLegendFontName',font);
            set(fig,'defaultTextFontName',font);
            set(fig, 'Units', 'inches');
            %set(fig, 'Position', [1, 1, 5, 3]); % Adjust figure size as per IEEE guidelines
            set(fig, 'Position', [1, 1, 5*1.25, 3*1.25]); % Adjust figure size as per IEEE guidelines
            set(ax, 'FontSize', 12); % Adjust font sizes as per IEEE guidelines
            set(ax, 'XGrid', 'on', 'YGrid', 'on');
            set(ax, 'YLim', [0, 0.1]);
            set(ax, 'XLim', [0, 100]);
            lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 12);
            set(lgd, 'Location', 'best');
            set(fig, 'Color', 'w');
            set(ax, 'Box','on');

            %title(sprintf('CQI = %d: Empirical vs. Fitted Gamma Distribution', cqi));
            xlabel('Downlink Bitrate (Mbps)');
            ylabel('Density');
            legend('Empirical Distribution', 'Fitted NB PDF');
            grid on;
            hold off;

            % export_fig NB_Fit_CQI_1.eps -m10
            % export_fig NB_Fit_CQI_1.png -m10
            % export_fig NB_Fit_CQI_1.pdf -m10
            % savefig( fig , 'NB_Fit_CQI_1.fig' )

            % Save figure as PDF
            saveas(fig, fullfile(output_dir, sprintf('NB_Fit_CQI_%d.pdf', cqi)));
            close(fig); % Close figure after saving to avoid memory issues
        else
            nbParams{cqi+1, 1} = cqi;
            nbParams{cqi+1, 2} = NaN; % Not enough variance for NB
            fprintf('CQI %2d: Variance too low, skipping NB fitting.\n', cqi);
        end
    else
        nbParams{cqi+1, 1} = cqi;
        nbParams{cqi+1, 2} = NaN; % Not enough data
        fprintf('CQI %2d: Not enough data to fit NB distribution.\n', cqi);
    end
end

% Save results to a .mat file
save('nb_params.mat', 'nbParams');
fprintf('\nNegative Binomial parameters saved to "nb_params.mat".\n');
