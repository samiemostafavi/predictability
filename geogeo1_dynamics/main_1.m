close all;
clear all;

% PLOT CLOSED FORM
K = 100;
Lres = 125;
Llims = [0,2000];
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
calc_and_plot_exact(3, c_mu, 0.85, K, Llims, Lres, zlims, "#545454", "-", "+");
hold on;
calc_and_plot_exact(3, c_mu, 0.75, K, Llims, Lres, zlims, "#545454", "-", "^");
hold on;

calc_and_plot_exact(9, c_mu, 0.85, K, Llims, Lres, zlims, "#545454", "-", "x");
hold on;
calc_and_plot_exact(9, c_mu, 0.75, K, Llims, Lres, zlims, "#545454", "-", "square");
hold on;

calc_and_plot_approx(9, c_mu, 0.85, K, Llims, Lres, zlims, "#174aa8", "--", "x");
hold on;
calc_and_plot_approx(9, c_mu, 0.75, K, Llims, Lres, zlims, "#174aa8", "--", "square");
hold on;

calc_and_plot_exact(15, c_mu, 0.85, K, Llims, Lres, zlims, "#545454", "-", "*");
hold on;
calc_and_plot_exact(15, c_mu, 0.75, K, Llims, Lres, zlims, "#545454", "-", "o");
hold on;

calc_and_plot_approx(15, c_mu, 0.85, K, Llims, Lres, zlims, "#174aa8", "--", "*");
hold on;
calc_and_plot_approx(15, c_mu, 0.75, K, Llims, Lres, zlims, "#174aa8", "--", "o");

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

% export_fig geogeo1_dynamics.eps -m10
% export_fig geogeo1_dynamics.png -m10
% savefig( fig , 'geogeo1_dynamics.fig' )


function result = calc_and_plot_exact(q_factor, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    expected_q = rho/(1-rho);
    q = ceil(q_factor*expected_q);
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, mm1exact(zlims, K, q, L, mu, lambda)];
    end
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Exact $x=%d \\chi, \\rho=%.2f$', q_factor, lambda/mu), 'LineWidth', 1.25);
    p.Color = color;
    p.LineStyle = linestyle;
    p.Marker = marker;
    p.MarkerSize = 10;
end

function result = calc_and_plot_approx(q_factor, c_mu, factor, K, Llims, Lres, zlims, color, linestyle, marker)
    res_cf = [];
    mu = c_mu;
    lambda = mu*factor;
    rho = lambda/mu;
    expected_q = rho/(1-rho);
    q = ceil(q_factor*expected_q);
    for L=[Llims(1):Lres:Llims(2)]
        res_cf = [res_cf, approx(q, L, mu, lambda)];
    end
    p = plot([Llims(1):Lres:Llims(2)],res_cf, '-o', 'DisplayName', sprintf('Approx $x=%d \\chi, \\rho=%.2f$', q_factor, lambda/mu), 'LineWidth', 1.25);
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