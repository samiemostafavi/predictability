close all;
clear all;

% PLOT CLOSED FORM
K = 128;
Lres = 100;
Llims = [0,1500];
zlims = [0,400];

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

%%%%%%%%%%%%%%%%%%%%%%%%%% PART ONE
c_mu = 0.5;
mu = c_mu;
util = 0.85;
lambda = c_mu*util;
rho = lambda/c_mu;
expected_q = rho/(1-rho)

P = sparse(K, K);
P(1, 1) = 1 - lambda * (1 - mu);
P(sub2ind([K, K], 2:K, 1:K-1)) = lambda * (1 - mu);
P(sub2ind([K, K], 1:K-1, 2:K)) = mu * (1 - lambda);
P(sub2ind([K, K], 2:K-1, 2:K-1)) = lambda * mu + (1 - lambda) * (1 - mu);
P(K, K) = 1 - mu * (1 - lambda);
P = full(P)';

tmp = P^10000;
pim = tmp(1,:);

% mm1 conditionals
conditionals = cell(K,1);
for n=[1:K]
    conditionals{n} = @(z) nbinpdf(z,n,mu);
end

fig = figure;
ax = axes(fig);
yplims = [16,64];
ypres = 16;
mylim = plot_conditionals(ax, conditionals,zlims,yplims,ypres);
plot_marginal_dist(ax, pim, conditionals, zlims);
comp_plot_conditionals(ax,fig,mylim,[0,100]); %35 for 0.8, 120 for 0.5

% export_fig geogeo1_validation2.eps -m10
% export_fig geogeo1_validation2.png -m10
% savefig( fig , 'geogeo1_validation2.fig' )

%%%%%%%%%%%%%%%%%%%%%%%%%% PART TWO

fig = figure;
ax = axes(fig);

calc_and_plot_exact(16, c_mu, util, K, Llims, Lres, zlims, "#545454", "-", "x", 1, 1);
hold on;
calc_and_plot_approx(16, c_mu, util, K, Llims, Lres, zlims, "#174aa8", "--", "x", 1, 1);
hold on;

calc_and_plot_exact(32, c_mu, util, K, Llims, Lres, zlims, "#545454", "-", "+", 1, 1);
hold on;
calc_and_plot_approx(32, c_mu, util, K, Llims, Lres, zlims, "#174aa8", "--", "+", 1, 1);
hold on;

calc_and_plot_exact(48, c_mu, util, K, Llims, Lres, zlims, "#545454", "-", "diamond", 1, 1);
hold on;
calc_and_plot_approx(48, c_mu, util, K, Llims, Lres, zlims, "#174aa8", "--", "diamond", 1, 1);
hold on;

calc_and_plot_exact(64, c_mu, util, K, Llims, Lres, zlims, "#545454", "-", "v", 1, 1);
hold on;
calc_and_plot_approx(64, c_mu, util, K, Llims, Lres, zlims, "#174aa8", "--", "v", 1, 1);
hold on;

%legend('closed form', 'simulation')

ylabel('Predictability');
xlabel('Lead Time');
font = 'Times New Roman';
set(fig,'defaultAxesFontName',font);
set(fig,'DefaultTextFontName', font, 'DefaultAxesFontName', font);
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

% export_fig geogeo1_evaluation2.eps -m10
% export_fig geogeo1_evaluation2.png -m10
% savefig( fig , 'geogeo1_evaluation2.fig' )


function comp_plot_conditionals(ax,fig,mylim,zlims)
    ylabel(ax, 'Probability');
    xlabel(ax, 'Sojourn Time Z');
    font = 'Times New Roman';
    set(fig,'defaultAxesFontName',font);
    set(fig,'DefaultTextFontName', font, 'DefaultAxesFontName', font);
    set(ax,'FontName', font, 'FontName', font);
    set(fig,'defaultLegendFontName',font);
    set(fig,'defaultTextFontName',font);
    set(fig, 'Units', 'inches');
    %set(fig, 'Position', [1, 1, 5, 3]); % Adjust figure size as per IEEE guidelines
    set(fig, 'Position', [1, 1, 5*1.25, 3*1.25/2]); % Adjust figure size as per IEEE guidelines
    set(ax, 'FontSize', 12); % Adjust font sizes as per IEEE guidelines
    set(ax, 'XGrid', 'on', 'YGrid', 'on');
    set(ax, 'YLim', [0, 0.25]);
    if length(zlims) ~= 0
        set(ax, 'XLim', [zlims(1), zlims(2)]);
    end
    lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 12);
    set(lgd, 'Location', 'best');
    set(fig, 'Color', 'w');
    set(ax, 'Box','on');
    %saveas(ax, 'final_conds_2.eps', 'epsc');
    %saveas(ax, 'final_conds_2.png');
end


function mylim = plot_conditionals(ax_cond, conditionals, zlims, yplims, ypres)
    
    legend_entry_added = false;
    % print conditionals
    K = length(conditionals);
    max_probs = [];
    for yp=[yplims(1):ypres:yplims(2)]
        if yp == 0
            y = 1;
        else
            y = yp;
        end
        probs = [];
        for z=[zlims(1):zlims(2)]
            probs = [probs, conditionals{y}(z)];
        end
        if sum(probs) < 0.99
           error('fix zlim for conditionals');
        end
        set(ax_cond, 'NextPlot', 'add');
        if ~legend_entry_added
            p = plot(ax_cond, [zlims(1):zlims(2)],probs, '--', 'LineWidth', 1, 'DisplayName', 'Analysis distributions');
            legend_entry_added = true;
        else
            p = plot(ax_cond, [zlims(1):zlims(2)],probs, '--', 'LineWidth', 1, 'HandleVisibility','off');
        end
        p.Color = "#545454";
        p.LineStyle = "--";
        
        max_probs = [max_probs, max(probs)];
    end
    
    % print conditionals
    K = length(conditionals);
    exp_vals = [];
    for y=[1:K]
        expected_value = 0;
        for z=[zlims(1):zlims(2)]
            expected_value = expected_value + z*conditionals{y}(z);
        end
        exp_vals = [exp_vals, expected_value];
    end
    
    i = 0;
    for yp=[yplims(1):ypres:yplims(2)]
        if yp == 0
            y = 1;
        else
            y = yp;
        end
        i = i+1;
        % Add text on top of each plot
        x_position = exp_vals(y) - 1.5;  % Middle of the x-axis range
        y_position = max_probs(i) * 1.1;          % Slightly above the maximum y value
        text(ax_cond, x_position, y_position, sprintf('x=%d', y), 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 12);
    end
    
    mylim = max_probs(1) * 1.1; 
    
    legend(ax_cond, 'show');
end

function plot_marginal_dist(ax, pim, conditionals, zlims)

    % print conditionals
    K = length(conditionals);
    probs = [];
    for z=[zlims(1):zlims(2)]
        prob = 0;
        for y=[1:K]
            prob = prob + pim(y)*conditionals{y}(z);
        end
        probs = [probs, prob];
    end
    p = plot(ax,[zlims(1):zlims(2)],probs, '-', 'LineWidth', 2.5, 'DisplayName', 'Marginal distribution');
    p.Color = "#545454";
    p.LineStyle = "-";
end

function result = calc_and_plot_exact(q, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker, offset, hop)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, mm1exact(zlims, K, q, L, mu, lambda)];
    end
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Exact $x=%d$', q), 'LineWidth', 1.25,'MarkerIndices',offset:hop:length(res_cf));
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function result = calc_and_plot_approx(q, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker, offset, hop)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, approx(q, L, mu, lambda)];
    end
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Approx. $x=%d$', q), 'LineWidth', 1.25,'MarkerIndices',offset:hop:length(res_cf));
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function result = approx(q, L, mu, lambda)

    beta = (lambda * (1 - mu)) / (mu * (1 - lambda));
    a = lambda * mu + (1 - lambda) * (1 - mu);
    b = 2 * sqrt(lambda * mu * (1 - lambda) * (1 - mu));
    c = b/(2*(a+b));

    result = (beta^((1-q)/2))*((1-sqrt(beta))^2)*((a+b)^L)/pi*theintegral(q, beta, L, c);
end

function result = theintegral(q, beta, L, c)  
    result = 0;
    dx = 0.0001;
    for x=[dx:dx:pi]
        result = result + sin(x)*sin(q*x)/((1-2*sqrt(beta)*cos(x)+beta)^2)*exp(-L*c*(x^2))*dx;
    end
end