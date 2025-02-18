close all;
clear all;

Llims = [1,1200]; % vehicle: 600, static: 1200
Lres = 20; % vehicle 10, static: 20
zlims = [0,1000];

%%%%%%%%%%%%%%%%%%%%%%%%%% PART ONE
%K = 64; % 64
%pr_stay = 0.2; % 0.2
%[P, pim] = get_lazy_random_walk_line_P(pr_stay,K);
%check_P(P, pim)
%conditionals = get_conditionals_linear_nbin(10,50,0.4,K);

K = 15; % 64
cqi_pr_stay = 0.90; % vehicle: 0.6 static: 0.9
[P, pim] = get_lazy_random_walk_line_P(cqi_pr_stay,K);
check_P(P, pim)
% Regression Model for Shape Parameter (k): k = 0.123 * CQI + 0.080
% Regression Model for Scale Parameter (θ): θ = 0.835 * CQI + 5.813
% alpha_k = 0.123;
% beta_k = 0.080;
% alpha_theta = 0.835;
% beta_theta = 5.813;
% conditionals = get_conditionals_gamma(alpha_k,beta_k,alpha_theta,beta_theta,K);

% Regression Model for r Parameter (r): r = 0.105 * CQI + 0.104
% Regression Model for p Parameter (p): p = -0.006 * CQI + 0.135
alpha_r = 0.105;
beta_r = 0.104;
alpha_p = -0.006;
beta_p = 0.136;
conditionals = get_conditionals_nbin(alpha_r,beta_r,alpha_p,beta_p,K);

% conds and marginals figure and its axes
fig = figure;
ax = axes(fig);
mylim = plot_conditionals(ax, conditionals,zlims);
plot_forecast_dist(ax, P, 60, 15, conditionals, zlims);
plot_marginal_dist(ax, pim, conditionals, zlims);
comp_plot_conditionals(ax,fig,mylim);


% export_fig randomwalk_validation_lte.eps -m10
% export_fig randomwalk_validation_lte.png -m10
% export_fig randomwalk_validation_lte.pdf -m10
% savefig( fig , 'randomwalk_validation_lte.fig' )


%%%%%%%%%%%%%%%%%%%%%%%%%% PART TWO

% results figure and its axes
fig = figure;
ax = axes(fig);

% colors: #3b00ff blue and #ed0105 red

x = K/15;
calc_and_plot_exact(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#545454", "-", "o", 1, 8);
calc_and_plot_upperbound1(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#174aa8", "--", "o", 1, 8);
calc_and_plot_upperbound2(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#e63538", ":", "o", 1, 8);

x = K/5;
calc_and_plot_exact(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#545454", "-", "^", 4, 8);
calc_and_plot_upperbound1(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#174aa8", "--", "^", 4, 8);
calc_and_plot_upperbound2(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#e63538", ":", "^", 4, 8);

x = K;
calc_and_plot_exact(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#545454", "-", "x", 6, 8);
calc_and_plot_upperbound1(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#174aa8", "--", "x", 6, 8);
calc_and_plot_upperbound2(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#e63538", ":", "x", 6, 8);

comp_plot_results(ax,fig,Llims);

% export_fig randomwalk_upperbounds_lte_vehicle.eps -m10
% export_fig randomwalk_upperbounds_lte_vehicle.png -m10
% export_fig randomwalk_upperbounds_lte_vehicle.pdf -m10
% savefig( fig , 'randomwalk_upperbounds_lte_vehicle.fig' )

% export_fig randomwalk_upperbounds_lte_static.eps -m10
% export_fig randomwalk_upperbounds_lte_static.png -m10
% export_fig randomwalk_upperbounds_lte_static.pdf -m10
% savefig( fig , 'randomwalk_upperbounds_lte_static.fig' )


function comp_plot_results(ax,fig,Llims)
    ylabel(ax, 'Predictability');
    xlabel(ax, 'Lead Time');
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
    set(ax, 'YLim', [0, 1]);
    set(ax, 'XLim', [Llims(1), Llims(2)]);
    lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 12);
    set(lgd, 'Location', 'best');
    set(fig, 'Color', 'w');
    set(ax, 'Box','on');
    %saveas(ax, 'final_main_2.eps', 'epsc');
    %saveas(ax, 'final_main_2.png');
end


function comp_plot_conditionals(ax,fig,mylim)
    ylabel(ax, 'Probability');
    xlabel(ax, 'Z [Throughput (Mbps)]');
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
    set(ax, 'YLim', [0, mylim]);
    set(ax, 'XLim', [0, 60]);
    lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 12);
    set(lgd, 'Location', 'best');
    set(fig, 'Color', 'w');
    set(ax, 'Box','on');
    %saveas(ax, 'final_conds_2.eps', 'epsc');
    %saveas(ax, 'final_conds_2.png');
end


function mylim = plot_conditionals(ax_cond, conditionals, zlims)
    
    legend_entry_added = false;
    HOP = 3;
    
    % Define marker styles
    markers = {'o', 's', '^', 'd', 'x', '+', '*', 'p', 'h'}; % Add more if needed
    num_markers = length(markers); % Total number of unique markers
    marker_idx = 1; % Initialize marker index
    marker_size = 10; % Set marker size
    offset =1;
    hop = 5;
    
    % print conditionals
    K = length(conditionals);
    max_probs = [];
    for yp=[0:HOP:K]
        if yp == 6 || yp == 9
           continue; 
        end
        if yp == 0
            y = 1;
        else
            y = yp;
        end
        
        % Extract probability values and z-values
        probs = [];
        z_values = [];
        for z = [zlims(1):zlims(2)]
            probs = [probs, conditionals{y}(z)];
            z_values = [z_values, z];
        end
        
        if sum(probs) < 0.99
           error('fix zlim for conditionals');
        end
        set(ax_cond, 'NextPlot', 'add');
        
        % Select a marker style (cycling through the list)
        marker_style = markers{mod(marker_idx - 1, num_markers) + 1};
    
        % Compute marker indices based on offset and hop
        marker_indices = offset:hop:length(z_values);

        % Plot the full curve
        p = plot(ax_cond, z_values, probs, '--', 'LineWidth', 1, ...
                 'DisplayName', sprintf('Posterior distribution for x=%d', y), ...
                 'Color', "#545454", 'Marker', marker_style, ...
                 'MarkerSize', marker_size, 'MarkerIndices', marker_indices);
        
        % Move to the next marker
        marker_idx = marker_idx + 1;
        
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
    for yp=[0:HOP:K]
        if yp == 0
            y = 1;
        else
            y = yp;
        end
        i = i+1;
        % Add text on top of each plot
        %x_position = exp_vals(y) - 1.5;  % Middle of the x-axis range
        %y_position = max_probs(i) * 1.05;          % Slightly above the maximum y value
        %text(ax_cond, x_position, y_position, sprintf('x=%d', y), 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 12);
    end
    
    %mylim = max_probs(1) * 1.1; 
    mylim = 0.08;
    
    legend(ax_cond, 'show');
end

function plot_forecast_dist(ax_cond, P, L, x, conditionals, zlims)
    lp = P^L;
    % print conditionals
    K = length(conditionals);
    probs = [];
    for z=[zlims(1):zlims(2)]
        prob = 0;
        for y=[1:K]
            prob = prob + lp(x,y)*conditionals{y}(z);
        end
        probs = [probs, prob];
    end
    p = plot(ax_cond, [zlims(1):zlims(2)],probs, '-.', 'LineWidth', 2.5, 'DisplayName', sprintf('Forecast distribution for L=%d, x=%d',L,x));
    p.Color = "#545454";
    p.LineStyle = "-.";
end

function plot_marginal_dist(ax_cond, pim, conditionals, zlims)
    
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
    p = plot(ax_cond, [zlims(1):zlims(2)],probs, '-', 'LineWidth', 2.5, 'DisplayName', 'Marginal distribution');
    p.Color = "#545454";
    p.LineStyle = "-";
end

function check_P(P, pim)

    % print pim
    fig_p = figure;
    ax_p = axes(fig_p);
    plot(ax_p,pim);

end

function conditionals = get_conditionals_nbin(alpha_r,beta_r,alpha_p,beta_p,K)

    % conditionals
    conditionals = cell(K,1);
    for y=[1:K]
        r = alpha_r*y + beta_r;
        p = alpha_p*y + beta_p;
        conditionals{y} = @(z) nbinpdf(z,r,p);
    end
end

function conditionals = get_conditionals_gamma(alpha_k,beta_k,alpha_theta,beta_theta,K)
    % conditionals
    % Regression Model for Shape Parameter (k): k = 0.123 * CQI + 0.080
    % Regression Model for Scale Parameter (θ): θ = 0.835 * CQI + 5.813
    conditionals = cell(K,1);
    for y=[1:K]
        k = alpha_k*y + beta_k;
        theta = alpha_theta*y + beta_theta;
        conditionals{y} = @(z) gampdf(z, k, theta);
    end
end

function [P,pim] = get_lazy_random_walk_line_P(p,K)
    % 1D walk on lattice
    q = (1 - p) / 2; % probability of moving left or right
    P = sparse(K, K);
    P(sub2ind([K, K], 1:K, 1:K)) = p;
    P(sub2ind([K, K], 2:K, 1:K-1)) = q; % from i to i+1
    P(sub2ind([K, K], 1:K-1, 2:K)) = q; % from i to i-1
    P(1, 1) = P(1, 1) + q; % stay at 1
    P(K, K) = P(K, K) + q; % stay at N
    P = full(P);
    
    tmp = P^100000;
    pim = tmp(1,:);
end

function calc_and_plot_exact(x, P, pim, conditionals, Llims, Lres, zlims, ax, color, linestyle, marker, offset, hop)
    K = length(P);
    results = [];
    for L=[Llims(1):Lres:Llims(2)]
        lp = P^L;
        l1norm = 0;
        for z=[zlims(1):zlims(2)]
            forecast = 0;
            marginal = 0;
            for y=[1:K]
                forecast = forecast + lp(x,y)*conditionals{y}(z);
                marginal = marginal + pim(y)*conditionals{y}(z);
            end
            l1norm = l1norm + abs(forecast-marginal);
        end
        result = 1/2*l1norm;
        results = [results, result];
    end
    set(ax, 'NextPlot', 'add');
    p = plot(ax,[Llims(1):Lres:Llims(2)],results, '-o', 'DisplayName', sprintf('Exact $x=%d$', x), 'LineWidth', 1.25, 'MarkerIndices',offset:hop:length(results));
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function calc_and_plot_upperbound2(x, P, pim, conditionals, Llims, Lres, zlims, ax, color, linestyle, marker, offset, hop)
    K = length(P);
    % Calc eigenvalues, take abs, and sort them
    [eigf, eigv] = weighted_eig(P,pim);
    eigv = diag(eigv);
    sorted_absed_eigs = sort(abs(eigv),'descend');
    eigstar = sorted_absed_eigs(2);
    
    % Calc R
    R=0;
    for z=[zlims(1):zlims(2)]
        y_arr = [];
        y2_arr = [];
        for y=[1:length(pim)]
            y_arr = [y_arr; conditionals{y}(z)];
            y2_arr = [y2_arr; (conditionals{y}(z)^2)];
        end
        R = R+((pim*y2_arr)/(pim*y_arr));
    end
    
    % Calc upperbound
    res_cf = [];   
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, (1/2)*(eigstar^L)*(((1/pim(x)-1))^0.5)*(sqrt(2)*((R-1)^0.5))];
    end
    set(ax, 'NextPlot', 'add');
    p = plot(ax,[Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB2 $x=%d$', x), 'LineWidth', 1.25, 'MarkerIndices',offset:hop:length(res_cf));
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function calc_and_plot_upperbound1(x, P, pim, conditionals, Llims, Lres, zlims, ax, color, linestyle, marker, offset, hop)
    K = length(P);
    % Calc eigenvalues, take abs, and sort them
    [eigf, eigv] = weighted_eig(P,pim);
    eigv = diag(eigv);
    sorted_absed_eigs = sort(abs(eigv),'descend');
    
    % Calc R
    R=0;
    for z=[zlims(1):zlims(2)]
        y_arr = [];
        y2_arr = [];
        for y=[1:length(pim)]
            y_arr = [y_arr; conditionals{y}(z)];
            y2_arr = [y2_arr; (conditionals{y}(z)^2)];
        end
        R = R+((pim*y2_arr)/(pim*y_arr));
    end
    
    % Calc upperbound
    res_cf = [];   
    K = length(pim);
    for L=[Llims(1):Lres:Llims(2)]
        mc_ub = 0;
        for i=[2:length(pim)]
            mc_ub = mc_ub + (eigf(x,i)^2)*(eigv(i)^(2*L));
        end
        res_cf = [res_cf, 0.5*(mc_ub^0.5)*(sqrt(2)*((R-1)^0.5))];
    end
    set(ax, 'NextPlot', 'add');
    p = plot(ax,[Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB1 $x=%d$', x), 'LineWidth', 1.25, 'MarkerIndices',offset:hop:length(res_cf));
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end