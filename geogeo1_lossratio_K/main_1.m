close all;
clear all;

% PLOT CLOSED FORM
Lres = 3;
Llims = [0,45];

% mu
c_mu = 0.35;

fig = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%

calc_and_plot_exact(1, c_mu, 0.9, 5, Llims, Lres, "#545454", "-", "o");
hold on;
calc_and_plot_exact(1, c_mu, 0.9, 9, Llims, Lres, "#545454", "-", "x");
hold on;
calc_and_plot_exact(1, c_mu, 0.9, 13, Llims, Lres, "#545454", "-", "square");
hold on;

calc_and_plot_exact(5, c_mu, 0.9, 5, Llims, Lres, "#545454", "-.", "o");
hold on;
calc_and_plot_exact(9, c_mu, 0.9, 9, Llims, Lres, "#545454", "-.", "x");
hold on;
calc_and_plot_exact(13, c_mu, 0.9, 13, Llims, Lres, "#545454", "-.", "square");
hold on;

ylabel('Predictability');
xlabel('Lead Time');
font = 'Times New Roman';
set(fig,'defaultAxesFontName',font);
set(fig,'DefaultTextFontName', font, 'DefaultAxesFontName', font);
ax = gca;
set(ax,'FontName', font, 'FontName', font);
set(fig,'defaultLegendFontName',font);
set(fig,'defaultTextFontName',font);
set(fig, 'Units', 'inches');
%set(fig, 'Position', [1, 1, 5, 3]); % Adjust figure size as per IEEE guidelines
set(fig, 'Position', [1, 1, 5*1.25, 3*1.25]); % Adjust figure size as per IEEE guidelines
set(gca, 'FontSize', 12); % Adjust font sizes as per IEEE guidelines
grid on;
legend('Interpreter', 'latex', 'FontSize', 12); % Use LaTeX for legend and adjust font size
legend('Location', 'best'); % Adjust legend location as per IEEE guidelines
set(gcf, 'Color', 'w');
ylim([0,1])




ax = gca;

% Get all line objects in current axes
hLines = findall(ax, 'Type', 'line');

% Define zoom area (example: around x=3 to x=5)
x1 = 30;
x2 = 45;
y1 = 0;
y2 = 0.05;

% Draw red rectangle on main plot
rectangle(ax, 'Position', [x1 y1 x2-x1 y2-y1], 'EdgeColor', 'r', 'LineWidth', 1.5);

% Create inset axes
axInset = axes('Position', [0.55 0.3 0.4 0.2]);
box on;
hold(axInset, 'on');

% Plot each line segment in the zoomed region
for i = 1:length(hLines)
    xData = get(hLines(i), 'XData');
    yData = get(hLines(i), 'YData');
    idx = xData >= x1 & xData <= x2;
    
    % Use original line's style
    plot(axInset, xData(idx), yData(idx), ...
        'Color', get(hLines(i), 'Color'), ...
        'LineStyle', get(hLines(i), 'LineStyle'), ...
        'LineWidth', get(hLines(i), 'LineWidth'), ...
        'Marker', get(hLines(i), 'Marker'), ...
        'MarkerSize', get(hLines(i), 'MarkerSize'), ...
        'MarkerEdgeColor', get(hLines(i), 'MarkerEdgeColor'), ...
        'MarkerFaceColor', get(hLines(i), 'MarkerFaceColor'));
end

xlim(axInset, [x1 x2]);
ylim(axInset, [y1 y2]);
set(axInset, 'XTick', [], 'YTick', []);

% Convert those to normalized figure coordinates
pt1 = axpos2figpos(ax, [x1 y1+y2]);  % top-left
pt2 = axpos2figpos(ax, [x2 y2]);  % top-right

% Define a point on the inset (you can refine these to your liking)
insetPos = get(axInset, 'Position');
insetTarget1 = [insetPos(1), insetPos(2)];  % bottom-left of inset
insetTarget2 = [insetPos(1) + insetPos(3), insetPos(2)];  % bottom-right of inset

% Add arrows or lines from corners of zoom box to inset
annotation('line', [pt1(1) insetTarget1(1)], [pt1(2) insetTarget1(2)], ...
    'Color', 'r', 'LineStyle', '--');
annotation('line', [pt2(1) insetTarget2(1)], [pt2(2) insetTarget2(2)], ...
    'Color', 'r', 'LineStyle', '--');


% export_fig geogeo1_lossratio_K.eps -m10
% export_fig geogeo1_lossratio_K.png -m10
% export_fig geogeo1_lossratio_K.pdf -m10
% savefig( fig , 'geogeo1_lossratio_K.fig' )


function result = calc_and_plot_exact(q, c_mu, factor, K, Llims, Lres, color, linestyle, marker)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    %expected_q = rho/(1-rho);
    expected_q = (rho/(1-rho)) - ((rho^(K+1) * (K+1)) / (1-(rho^(K+1))));
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, mm1exactlossratio(K, q, L, mu, lambda)];
    end
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Exact $x=%d, \\rho=%.2f, K=%d$', q, lambda/mu, K), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end
