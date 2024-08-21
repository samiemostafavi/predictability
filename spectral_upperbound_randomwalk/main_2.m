close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
Llims = [1,3000];
Lres = 150;
zlims = [0,150];
fig = figure;
ax = axes(fig,'Position', [0.1 0.6 0.8 0.35]);

K = 64;
pr_stay = 0.1;
[P, pim] = get_lazy_random_walk_line_P(pr_stay,K); 
check_P(P, pim)
conditionals = get_conditionals_linear_nbin(10,50,0.4,K);
check_conditionals(conditionals,zlims)
x = K/2;
calc_and_plot_exact(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#545454", "-", "s", {'x', 'K'});
calc_and_plot_upperbound1(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#174aa8", "--", "s", {'x', 'K'});
calc_and_plot_upperbound2(x, P, pim, conditionals, Llims, Lres, zlims, ax, "#e63538", ":", "s", {'x', 'K'});
doplot(fig,ax)

%%%%%%%%%%%%%%%%%%%%%%%%%%
Llims = [1,810];
Lres = 38;
zlims = [0,150];

ax = axes(fig,'Position', [0.1 0.15 0.8 0.35]);

Kp = 16;
[aggP, aggConditionals] = aggregateMMP(P, pim, conditionals, Kp);
% form stationary matrix
tmp = aggP^10000;
aggPim = tmp(1,:);
check_P(aggP, aggPim)
check_conditionals(aggConditionals,zlims)
x = Kp/2;
calc_and_plot_exact(x, aggP, aggPim, aggConditionals, Llims, Lres, zlims, ax, "#545454", "-", "*", {'a', '\bar{K}'});
calc_and_plot_upperbound1(x, aggP, aggPim, aggConditionals, Llims, Lres, zlims, ax, "#174aa8", "--", "*", {'a', '\bar{K}'});
calc_and_plot_upperbound2(x, aggP, aggPim, aggConditionals, Llims, Lres, zlims, ax, "#e63538", ":", "*", {'a', '\bar{K}'});
doplot(fig,ax)
set(ax, 'XTick', [0:135:Llims(2)+1])
xlabel(ax, 'Lead Time');
ylab = ylabel(ax, 'Predictability');
ylab.Position(2) = 0.6; % Vertically center it in the figure
 
% export_fig randomwalk_multiple_upperbounds.eps -m10
% export_fig randomwalk_multiple_upperbounds.png -m10
% savefig( fig , 'randomwalk_multiple_upperbounds.fig' )

function doplot(fig,ax)

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
    set(ax, 'YLim', [0, 0.5]);
    lgd = legend(ax, 'Interpreter', 'latex', 'FontSize', 12);
    set(lgd, 'Location', 'best');
    set(fig, 'Color', 'w');
    set(ax, 'Box','on');
end

function check_conditionals(conditionals, zlims)

    % print conditionals
    fig_cond = figure;
    ax_cond = axes(fig_cond);
    K = length(conditionals);
    for y=[1:K]
        probs = [];
        for z=[zlims(1):zlims(2)]
            probs = [probs, conditionals{y}(z)];
        end
        if sum(probs) < 0.99
           error('fix zlim for conditionals');
        end
        hold on;
        plot(ax_cond,[zlims(1):zlims(2)],probs);
    end
    
    
    % print conditionals
    fig_cond = figure;
    ax_cond = axes(fig_cond);
    K = length(conditionals);
    exp_vals = [];
    for y=[1:K]
        expected_value = 0;
        for z=[zlims(1):zlims(2)]
            expected_value = expected_value + z*conditionals{y}(z);
        end
        exp_vals = [exp_vals, expected_value];
    end
    plot(ax_cond,[1:K], exp_vals, '-*');
    xlim([1,K])
end

function check_P(P, pim)

    % print pim
    fig_p = figure;
    ax_p = axes(fig_p);
    plot(ax_p,pim);

end

function conditionals = get_conditionals_linear_nbin(r1, r2, p, K)
    % conditionals
    conditionals = cell(K,1);
    for y=[1:K]
        r = r1 + (y-1)*(r2-r1)/K;
        conditionals{y} = @(z) nbinpdf(z,r,p);
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
    P = full(P)';
    
    tmp = P^100000;
    pim = tmp(1,:);
end

function calc_and_plot_exact(x, P, pim, conditionals, Llims, Lres, zlims, ax, color, linestyle, marker, legends)
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
    p = plot(ax,[Llims(1):Lres:Llims(2)],results, '-o', 'DisplayName', sprintf('Exact $%s=%d, %s=%d$', legends{1}, x, legends{2}, K), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function calc_and_plot_upperbound2(x, P, pim, conditionals, Llims, Lres, zlims, ax, color, linestyle, marker, legends)
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
        res_cf = [res_cf, (1/2)*(eigstar^L)*(((1/pim(x)-1))^0.5)];
    end
    set(ax, 'NextPlot', 'add');
    p = plot(ax,[Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB2 $%s=%d, %s=%d$', legends{1}, x, legends{2}, K), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function calc_and_plot_upperbound1(x, P, pim, conditionals, Llims, Lres, zlims, ax, color, linestyle, marker, legends)
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
        res_cf = [res_cf, 0.5*(mc_ub^0.5)*((R-1)^0.5)];
    end
    set(ax, 'NextPlot', 'add');
    p = plot(ax,[Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB1 $%s=%d, %s=%d$', legends{1}, x, legends{2}, K), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end