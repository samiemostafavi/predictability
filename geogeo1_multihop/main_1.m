close all;
clear all;

% PLOT CLOSED FORM
K = 100;
Llims = [0,2000];
Lres = 150;
zlims = [0,600];

% hops
M = 5;

% mu
c_mus = [0.4, 0.4, 0.4, 0.4, 0.4];

% multihop mm1 conditionals
conditionals = cell(K,M);
for m=[1:M]
    for n=[1:K]
        conditionals{n,m} = @(z) nbinpdf(z,n,c_mus(m));
    end
end

% print conditionals
% tiledlayout(M,1);
% K = length(conditionals);
% for m=[1:M]
%     nexttile
%     for n=[1:K]
%         probs = [];
%         for z=[zlims(1):zlims(2)]
%             probs = [probs, conditionals{n,m}(z)];
%         end
%         hold on;
%         plot([zlims(1):zlims(2)],probs);
%         title(sprintf('Hop %d',m));
%     end
% end

% fig = figure;
fig = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%

Ks = [100, 100, 100, 100, 100];
first_util = 0.8;


q_factors = [15,15,15,15,15];
calc_and_plot_multi(q_factors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#545454", "-", "o");

hold on;
q_factor = 15;
K = Ks(1);
c_mu = c_mus(1);
calc_and_plot_single_times(5, q_factor, c_mu, first_util, K, Llims, Lres, zlims, "#174aa8", "--", "o");

hold on;
q_factors = [15,15,15,-1,-1];
calc_and_plot_multi(q_factors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#545454", "-", "x");

hold on;
q_factors = [15,-1,-1,-1,-1];
calc_and_plot_multi(q_factors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#545454", "-", "^");

hold on;
q_factor = 15;
K = Ks(1);
c_mu = c_mus(1);
calc_and_plot_single(q_factor, c_mu, first_util, K, Llims, Lres, zlims, "#174aa8", "--", "^",1);


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

% export_fig geogeo1_multihop1.eps -m10
% export_fig geogeo1_multihop1.png -m10
% savefig( fig , 'geogeo1_multihop1.fig' )

function result = calc_and_plot_single_times(times, q_factor, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    expected_q = rho/(1-rho);
    q = ceil(q_factor*expected_q);
    for L=[Llims(1):Lres:Llims(2)]
        % mm1closedform(ylims, yres, K, q, L, mu, lambda)
        res_cf = [res_cf, times*mm1exact(zlims, K, q, L, mu, lambda)];
    end
    %plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('$q=%d, \\mu=%.2f, \\lambda=%.2f$', q,mu,lambda), 'LineWidth', 1)
    %p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('$\\sum_{m=1}^{5} S^{(m)}_{L,%d\\chi^{(m)}}$',q_factor), 'LineWidth', 1);
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB $x=[%s]$',array2strf2([q_factor,q_factor,q_factor,q_factor,q_factor])), 'LineWidth', 1.25);

    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function result = calc_and_plot_single(q_factor, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker, m)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    expected_q = rho/(1-rho);
    q = ceil(q_factor*expected_q);
    for L=[Llims(1):Lres:Llims(2)]
        % mm1closedform(ylims, yres, K, q, L, mu, lambda)
        res_cf = [res_cf, mm1exact(zlims, K, q, L, mu, lambda)];
    end
    %plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('$q=%d, \\mu=%.2f, \\lambda=%.2f$', q,mu,lambda), 'LineWidth', 1)
    %p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('$x^{(%d)}_n=%d\\chi^{(%d)}$', m, q_factor, m), 'LineWidth', 1);
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB $x=[%s]$',array2strf2([q_factor,-1,-1,-1,-1])), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function result = calc_and_plot_multi(q_factors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, color, linestyle, marker)
    
    pims = [];
    Ps = cell(M);
    lambda = c_mus(1)*first_util;
    qs = [];
    
    for idx=[1:M]
        mu = c_mus(idx);
        K = Ks(idx);
        
        rho = lambda/mu;
        expected_q = rho/(1-rho);
        if q_factors(idx) > 0
            q = ceil(q_factors(idx)*expected_q);
            qs = [qs, q];
        else
            qs = [qs, -1];
        end
            
        % form m/m/1 P matrix
        P = sparse(K, K);
        P(1, 1) = 1 - lambda * (1 - mu);
        P(sub2ind([K, K], 2:K, 1:K-1)) = lambda * (1 - mu);
        P(sub2ind([K, K], 1:K-1, 2:K)) = mu * (1 - lambda);
        P(sub2ind([K, K], 2:K-1, 2:K-1)) = lambda * mu + (1 - lambda) * (1 - mu);
        P(K, K) = 1 - mu * (1 - lambda);
        P = full(P)';
        Ps{idx} = P;
        
        tmp = P^10000;
        pim = tmp(1,:);
        pims = [pims; pim];
    end
    
    res_cf = [];
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, predictability_exact_multi(zlims, qs, L, Ps, pims, M, Ks, conditionals)];
    end
    % p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('$x_n=[%s], \\mu=[%s], \\lambda=%.2f$',array2str(qs), array2strf(c_mus), lambda), 'LineWidth', 1);
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Exact $x=[%s]$',array2strf2(q_factors)), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function str = array2strf(vec)
    str = sprintf(',%.2f',vec);
    str = str(2:end);
end

function str = array2strf2(vec)
    str = '';
    for i=[1:length(vec)]
        ch = '';
        if vec(i) < 0
           ch = '-';
        else
           ch = sprintf('%d\\chi',vec(i));
        end
        str = strcat(str,ch,',');
    end
    str = str(1:end-1);
end

function str = array2str(vec)
    str = '';
    for i=[1:length(vec)]
        ch = '';
        if vec(i) < 0
           ch = '-'; 
        else
           ch = sprintf('%d',vec(i));
        end
        str = strcat(str,ch,',');
    end
    str = str(1:end-1);
end