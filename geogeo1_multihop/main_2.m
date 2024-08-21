close all;
clear all;

% PLOT CLOSED FORM
K = 100;
Llims = [0,2000];
Lres = 80;
zlims = [0,600];

% hops
M = 3;

% mu
c_mus = [0.4, 0.38, 0.4];

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

fig = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%

Ks = [100, 100, 100];
first_util = 0.8;
mid_util = first_util*c_mus(1)/c_mus(2);
expected_first_q = first_util/(1-first_util);
expected_mid_q = mid_util/(1-mid_util);

qfactors = [10,10,10];
qs = [ceil(expected_first_q*qfactors(1)),ceil(expected_mid_q*qfactors(2)),ceil(expected_first_q*qfactors(3))];
calc_and_plot(qs, qfactors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#545454", "-", "^", 1, 3);

hold on;
qfactors = [10,10,10];
qs = [ceil(expected_first_q*qfactors(1)),ceil(expected_mid_q*qfactors(2)),ceil(expected_first_q*qfactors(3))];
calc_and_plot_add_multi(qs, qfactors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#174aa8", "--", "^", 1, 3);

hold on;
qs = [ceil(expected_first_q*10),-1,-1];
qfactors = [10,-1,-1];
calc_and_plot(qs, qfactors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#545454", "-", "o", 1, 3);

hold on;
qs = [-1,ceil(expected_mid_q*10),-1];
qfactors = [-1,10,-1];
calc_and_plot(qs, qfactors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#545454", "-", "s", 2, 3);

hold on;
qs = [-1,-1,ceil(expected_first_q*10)];
qfactors = [-1,-1,10];
calc_and_plot(qs, qfactors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, "#545454", "-", "x", 3, 3);

hold on;
q_factor = 10;
K = Ks(1);
c_mu = c_mus(1);
qfactors = [10,-1,-1];
calc_and_plot_single(q_factor, qfactors, c_mu, first_util, K, Llims, Lres, zlims, "#174aa8", "--", "o", 1, 3);

hold on;
q_factor = 10;
K = Ks(2);
c_mu = c_mus(2);
qfactors = [-1,10,-1];
calc_and_plot_single(q_factor, qfactors, c_mu, mid_util, K, Llims, Lres, zlims, "#174aa8", "--", "s", 1, 3);


% export_fig geogeo1_multihop2.eps -m10
% export_fig geogeo1_multihop2.png -m10
% savefig( fig , 'geogeo1_multihop2.fig' )

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

saveas(gca, 'final_closedformmulti_2.eps', 'epsc');
saveas(gca, 'final_closedformmulti_2.png');

% export_fig final_closedformmulti_2.png -m10


function result = calc_and_plot(qs, qfactors, c_mus, first_util, Ks, M, conditionals, Llims, Lres, zlims, color, linestyle, marker, offset, hop)
    
    pims = [];
    Ps = cell(M);
    lambda = c_mus(1)*first_util;
    
    for idx=[1:M]
        q = qs(idx);
        mu = c_mus(idx);
        K = Ks(idx);
       
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
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Exact $x=[%s]$',array2str(qfactors)), 'LineWidth', 1.25, 'MarkerIndices', offset:hop:length(res_cf));
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function result = calc_and_plot_add_multi(qs, qfactors, c_mus, factor, Ks, M, conditionals, Llims, Lres, zlims, color, linestyle, marker, offset, hop)
    res_cf = [];
    lambda = c_mus(1)*factor;

    for L=[Llims(1):Lres:Llims(2)]
        qres = 0;
        for qi=[1:M]
            mu = c_mus(qi);
            q = qs(qi);
            K = Ks(qi);
            qres = qres + mm1exact(zlims, K, q, L, mu, lambda);
        end
        res_cf = [res_cf, qres];
    end
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB $x=[%s]$',array2str(qfactors)), 'LineWidth', 1.25, 'MarkerIndices', offset:hop:length(res_cf));

    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function result = calc_and_plot_single(q_factor, qfactors, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker, offset, hop)
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
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('UB $x=[%s]$',array2str(qfactors)), 'LineWidth', 1.25, 'MarkerIndices', offset:hop:length(res_cf));
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function str = array2strf(vec)
    str = sprintf(',%.2f',vec);
    str = str(2:end);
end

function str = array2str(vec)
    str = '';
    for i=[1:length(vec)]
        ch = '';
        if vec(i) < 0
           ch = '-'; 
        else
           ch = sprintf('%d \\chi^{(%d)}',vec(i),i);
        end
        str = strcat(str,ch,',');
    end
    str = str(1:end-1);
end