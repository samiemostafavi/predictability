close all;
clear all;

% PLOT CLOSED FORM
Lres = 40;
Llims = [0,400];
zlims = [0,400];

% mu
c_mu = 0.4;

% print conditionals
% figure;
% K = length(conditionals);
% for i=[1:K]
%     probs = [];
%     for z=[zlims(1):zlims(2)]
%         probs = [probs, nbinpdf(z,conditionals(i,1),conditionals(i,2))];
%     end
%     hold on;
%     plot([zlims(1):zlims(2)],probs);
% end

fig = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%

calc_and_plot_exact(3, c_mu, 0.85, 12, Llims, Lres, zlims, "#545454", "-", "o");
hold on;
calc_and_plot_exact(3, c_mu, 0.85, 15, Llims, Lres, zlims, "#545454", "-", "^");
hold on;
calc_and_plot_exact(3, c_mu, 0.85, 20, Llims, Lres, zlims, "#545454", "-", "*");
hold on;
calc_and_plot_exact(3, c_mu, 0.85, 25, Llims, Lres, zlims, "#545454", "-", "x");
hold on;
calc_and_plot_exact(3, c_mu, 0.85, 30, Llims, Lres, zlims, "#545454", "-", "square");
hold on;
calc_and_plot_exact(3, c_mu, 0.85, 128, Llims, Lres, zlims, "#545454", "-", "+");
hold on;

%legend('closed form', 'simulation')

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

% export_fig geogeo1_dynamics_K.eps -m10
% export_fig geogeo1_dynamics_K.png -m10
% export_fig geogeo1_dynamics_K.pdf -m10
% savefig( fig , 'geogeo1_dynamics_K.fig' )


function result = calc_and_plot_exact(q_factor, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    %expected_q = rho/(1-rho);
    expected_q = (rho/(1-rho)) - ((rho^(K+1) * (K+1)) / (1-(rho^(K+1))));
    q = ceil(q_factor*expected_q);
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, mm1exact(zlims, K, q, L, mu, lambda)];
    end
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Exact $x=%d \\chi, \\rho=%.2f, K=%d$', q_factor, lambda/mu, K), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end
