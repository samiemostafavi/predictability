close all;
clear all;

% PLOT CLOSED FORM
K = 50;
Llims = [1,400];
Lres = 5;
qlims = [1,50];
zlims = [0,300]; %250

% mu
c_mu = 0.2;

% print conditionals
% conditionals = [];
% for n=[1:K]
%     conditionals = [conditionals; [n, c_mu]];
% end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%

% test mMUs
mulims = [0.2,0.5];
mures = 0.05;
mus = [mulims(1):mures:mulims(2)];

% Notice: running each of the following blocks takes 1 hr, the results are saved in files

% util = 0.925;
% resp925 = [];
% for mu = [mulims(1):mures:mulims(2)]
%     resp925 = [resp925, calc_predtime(0.1, mu, util, K, Llims, Lres, zlims, qlims)];
% end
% resp925 = [resp925;mus];

% util = 0.9;
% resp9 = [];
% for mu = [mulims(1):mures:mulims(2)]
%     resp9 = [resp9, calc_predtime(0.1, mu, util, K, Llims, Lres, zlims, qlims)];
% end
% resp9 = [resp9;mus];

% util = 0.875;
% resp875 = [];
% for mu = [mulims(1):mures:mulims(2)]
%     resp875 = [resp875, calc_predtime(0.1, mu, util, K, Llims, Lres, zlims, qlims)];
% end
% resp875 = [resp875;mus];

% util = 0.85;
% resp85 = [];
% for mu = [mulims(1):mures:mulims(2)]
%     resp85 = [resp85, calc_predtime(0.1, mu, util, K, Llims, Lres, zlims, qlims)];
% end
% resp85 = [resp85;mus];

% load from files
load('predtime_utilp875.mat')
load('predtime_utilp85.mat')
load('predtime_utilp9.mat')
load('predtime_utilp925.mat')

fig = figure;

p = plot(mus,resp85(1,:), '-o', 'DisplayName', sprintf('$\\rho=%.3f$', 0.85), 'LineWidth', 1.25);
p.Color = "#545454";
p.LineStyle = "-";
p.Marker = "o";
p.MarkerSize = 10;
hold on;

p = plot(mus,resp875(1,:), '-o', 'DisplayName', sprintf('$\\rho=%.3f$', 0.875), 'LineWidth', 1.25);
p.Color = "#545454";
p.LineStyle = "-";
p.Marker = "^";
p.MarkerSize = 10;
hold on;

p = plot(mus,resp9(1,:), '-o', 'DisplayName', sprintf('$\\rho=%.3f$', 0.9), 'LineWidth', 1.25);
p.Color = "#545454";
p.LineStyle = "-";
p.Marker = "*";
p.MarkerSize = 10;

p = plot(mus,resp925(1,:), '-o', 'DisplayName', sprintf('$\\rho=%.3f$', 0.925), 'LineWidth', 1.25);
p.Color = "#545454";
p.LineStyle = "-";
p.Marker = "s";
p.MarkerSize = 10;

ylabel('Epsilon-predictable horizon');
xlabel('Service Rate');
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

% export_fig geogeo1_predtime.eps -m10
% export_fig geogeo1_predtime.png -m10
% savefig( fig , 'geogeo1_predtime.fig' )

function result = calc_predtime(epsilon, mu, factor, K, Llims, Lres, zlims, qlims)
    lambda = mu*factor;
    rho = lambda/mu;
    expected_q = rho/(1-rho);
    %q = ceil(q_factor*expected_q);
    res_cf = [];
    for L=[Llims(1):Lres:Llims(2)]
        res_q = [];
        for q=[qlims(1):qlims(2)]
            res_q = [res_q, mm1exact(zlims, K, q, L, mu, lambda)];
        end
        res_cf = [res_cf, min(res_q)];
    end
    res_idx = find(diff(res_cf<=epsilon));
    Ls = [Llims(1):Lres:Llims(2)];
    result = Ls(res_idx);
end