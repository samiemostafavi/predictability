clc; clear; close all;

% Load the previously saved Gamma parameters
load('nb_params.mat', 'nbParams');

% Extract CQI values and Gamma parameters
cqi_values = cell2mat(nbParams(:,1)); % CQI levels

% Ensure gamma parameters are extracted correctly
validNB = ~cellfun(@isempty, nbParams(:,2)) & ~cellfun(@(x) any(isnan(x)), nbParams(:,2));
validNBParams = nbParams(validNB, 2); % Keep only valid rows

% Convert to matrix
nb_params = cell2mat(validNBParams); % Shape (k) and Scale (θ)

% Filter CQI values accordingly
cqi_values = cqi_values(validNB);

% Separate shape (k) and scale (θ)
r_values = nb_params(:,1);
p_values = nb_params(:,2);

% Perform linear regression
r_fit = polyfit(cqi_values, r_values, 1); % Linear fit for (r)
p_fit = polyfit(cqi_values, p_values, 1); % Linear fit for (p)

% Generate fitted values
fitted_r = polyval(r_fit, cqi_values);
fitted_p = polyval(p_fit, cqi_values);

% Display regression results
fprintf("Regression Model for r Parameter (r): r = %.3f * CQI + %.3f\n", r_fit(1), r_fit(2));
fprintf("Regression Model for p Parameter (p): p = %.3f * CQI + %.3f\n", p_fit(1), p_fit(2));

% Plot Shape (k) vs CQI

% Create plot
fig = figure();
ax = axes(fig);

scatter(ax, cqi_values, r_values, 'b', 'filled'); hold on;
plot(ax, cqi_values, fitted_r, 'r-', 'LineWidth', 2);

%title('Linear Regression: r Parameter vs CQI');
%legend('Observed', 'Fitted', 'Location', 'NorthWest');
%grid on;
%hold off;
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
lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 12);
set(lgd, 'Location', 'best');
set(fig, 'Color', 'w');
set(ax, 'Box','on');

xlabel('CQI');
ylabel('NB r Parameter (r)');
legend('Observed', 'Fitted', 'Location', 'NorthWest');
grid on;
hold off;

% export_fig NB_regression_CQI_r.eps -m10
% export_fig NB_regression_CQI_r.png -m10
% export_fig NB_regression_CQI_r.pdf -m10
% savefig( fig , 'NB_regression_CQI_r.fig' )


% Plot Scale (θ) vs CQI

% Create plot
fig = figure();
ax = axes(fig);
scatter(ax, cqi_values, p_values, 'g', 'filled'); hold on;
plot(ax, cqi_values, fitted_p, 'r-', 'LineWidth', 2);
xlabel('CQI');
ylabel('NB p Parameter (p)');

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
lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 12);
set(lgd, 'Location', 'best');
set(fig, 'Color', 'w');
set(ax, 'Box','on');

%title('Linear Regression: p Parameter vs CQI');
legend('Observed', 'Fitted', 'Location', 'NorthWest');
grid on;
hold off;

% export_fig NB_regression_CQI_p.eps -m10
% export_fig NB_regression_CQI_p.png -m10
% export_fig NB_regression_CQI_p.pdf -m10
% savefig( fig , 'NB_regression_CQI_p.fig' )


% Save results
r_p_fit = [r_fit; p_fit];
save('cqi_bitrate_regression_nb.mat', 'r_p_fit');
fprintf('\nRegression results saved to "cqi_bitrate_regression_nb.mat".\n');
