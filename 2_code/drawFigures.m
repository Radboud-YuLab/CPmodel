function [Fig1d, Fig2a, Fig2b, Fig3, Fig4a, Fig4b] = drawFigures 
% Draw all figures in the CP paper.
% required input: pom_main.mat, uhl_main.mat, and uhl_pred_CYHX.mat (from 
% the script CP_paper_main.m.
% output: 
%       Fig1d - best fit models a single gene (Il4i1 from
%       the Pommerenke dataset)
%       Fig2a - Pommenrenke 2012 CP density plots
%       Fig2b - Uhlitz 2017 CP density plots
%       Fig3 - predictions of the Uhlitz CYHX dataset
%       Fig4a - Qu 2018 CP density plots, 
%       Fig4b - predictions of EEC syndrome in Qu 2018

% Chen chen. Last update: 2024-09-16
% Rosemary Yu. Last update: 2025-02-08


%% Figure 1d - best fit model of a single gene from the Pommerenke dataset
load('4_processed_data\pom_main\pom_main.mat');
idx = find (strcmp ('Il4i1', out_opt (:,1)));  
output = out_all(idx);
nCPs = output.n_CPs;
coords = output.CP_coordinates;
data = pom_log2_fc (idx,:);
data = transpose (cat (1, pom_days, data));

Fig1d = figure;
hold on
p1 = scatter (data (:,1), data (:,2), 10,"black","filled");

P = graph_points (nCPs (1,:), coords {1,:});
p2 = plot(P(:,1), P(:,2),"Color",[0 0.64 0.93],"LineWidth",2.0);

P = graph_points (nCPs (2,:), coords {2,:});
p3 = plot(P(:,1), P(:,2),"Color",[0.8500 0.3250 0.0980],"LineWidth",2.0);

P = graph_points (nCPs (3,:), coords {3,:});
p4 = plot(P(:,1), P(:,2),"Color",[0.4660 0.6740 0.1880],"LineWidth",2.0);

P = graph_points (nCPs (4,:), coords {4,:});
p5 = plot(P(:,1), P(:,2),"Color",[0.9290 0.6940 0.1250],"LineWidth",2.0);

xlabel ('time (days)')
ylabel ('gene expression (log2FC)')
ylim([0 3])
legend ([p1 p2 p3 p4 p5], ...
    {'data', 'model, n = 1', 'model, n = 2', 'model, n = 3', 'model, n = 4'});

RSS = output.RSS;
nCP_penalized_RSS = output.nCP_penalized_RSS;
R2 = output.R2;
save('4_processed_data\pom_main\Fig1d_table.mat', 'RSS', 'nCP_penalized_RSS', 'R2');
saveas(Fig1d,'6_results\figures\Figure1d.png')

%% Figure 2a - Pommenrenke 2012 CP density plots
Fig2a = figure;
hold on

%B cell response
[cluster1_2] = findCluster ([1;2], out_opt, pom_clusters);
[f,xi] = drawDensity(cluster1_2, pom_days, 0.1);
f = f.*0.18; 
area(xi,f, "FaceColor", [0 0.64 0.93], "FaceAlpha", 0.2, "EdgeColor", "none");
d2 = plot(xi,f, "Color", [0 0.64 0.93], "LineWidth", 1.5);

% T cell response
[cluster4] = findCluster (4, out_opt, pom_clusters);
[f,xi] = drawDensity(cluster4, pom_days, 0.1);
f = f.*0.45; 
area(xi,f, "FaceColor", [0.8500 0.3250 0.0980], "FaceAlpha", 0.2, "EdgeColor", "none");
d3 = plot(xi,f, "Color",[0.8500 0.3250 0.0980], "LineWidth", 1.5);

% innate
[cluster6] = findCluster (6, out_opt, pom_clusters);
[f,xi] = drawDensity(cluster6, pom_days, 0.1);
f = f.*0.17; 
area(xi,f, "FaceColor", [0.4660 0.6740 0.1880], "FaceAlpha", 0.2, "EdgeColor", "none");
d4 = plot(xi,f, "Color",[0.4660 0.6740 0.1880], "LineWidth", 1.5);

% repair
[cluster7] = findCluster (7, out_opt, pom_clusters);
[f,xi] = drawDensity(cluster7, pom_days, 0.1);
f = f.*0.3; 
area(xi,f, "FaceColor", [0.9290 0.6940 0.1250], "FaceAlpha", 0.2, "EdgeColor", "none");
d5 = plot(xi,f, "Color",[0.9290 0.6940 0.1250], "LineWidth", 1.5);

%all data
[f,xi] = drawDensity(out_opt, pom_days, 0.1);
d1 = plot(xi,f,"Color",'black',"LineWidth", 2.0);

daspect([25 0.1 1]) % 0.5*xlim, 1*ylim, and 1
xlabel ('time (days post infection)');
ylabel ('density');
legend ([d1 d4 d3 d2 d5], {'all data', 'innate response', 'T cell response', 'B cell response', 'tissue repair'});

saveas(Fig2a,'6_results\figures\Figure2a.png')


%% Figure 2b - Uhlitz 2017 CP density plots
load('4_processed_data\uhl_main\uhl_main.mat');
Fig2b = figure;
hold on

[IEG_cluster] = findCluster ('IEG', out_opt, uhl_cluster);
[f,xi] = drawDensity(IEG_cluster, uhl_hrs, 0.2);
f = f.*0.055; %*0.015;
area(xi,f, "FaceColor", [0.85 0.325 0.098], "FaceAlpha", 0.2, "EdgeColor", "none");
g2 = plot(xi,f, "Color",[0.85 0.325 0.098],"LineWidth",1.5);

[SRG_cluster] = findCluster ('SRG', out_opt, uhl_cluster);
[f,xi] = drawDensity(SRG_cluster, uhl_hrs, 0.2);
f = f.*0.7; %0.45;
area(xi,f, "FaceColor", [0 0.64 0.93],"FaceAlpha", 0.25, "EdgeColor", "none");
g3 = plot(xi,f, "Color",[0 0.64 0.93],"LineWidth",1.5);

[f,xi] = drawDensity(out_opt, uhl_hrs, 0.2);
g1 = plot(xi,f, "Color", 'black',"LineWidth",2.0);

tf = {'EGR1'; 'FOS'; 'FOSB'; 'EGR2'; 'JUNB'; 'EGR3'; 'KLF10'; 'ARID3B';...
    'TFB2M'; 'FOSL1'; 'ZFP36'; 'FOSL2'; 'EGR4'; 'NR4A1'};
tf_xi = [];
for i = 1:height(tf)
    idx = find(ismember(out_opt(:,1), tf(i)));
    cp1 = cell2mat(out_opt(idx,3));
    tf_xi = [tf_xi; cp1(2,1)];
end
tf_yi(1:14,1)=0.5;
gtf1 = plot(tf_xi(1:6), tf_yi(1:6), '|', 'Color',[0.85 0.325 0.098],...
    "MarkerSize", 10, "LineWidth", 1.2);
gtf2 = plot(tf_xi(7:14), tf_yi(7:14), '|', 'Color',[0 0.64 0.93],...
    "MarkerSize", 10, "LineWidth", 1.2);
gtf = plot(0,-1,'|', 'Color', 'black', "MarkerSize", 10, "LineWidth", 1.2);

daspect([4 0.6 1]) % 0.5*xlim, 1*ylim, and 1
ylim([0 0.6])
xlabel ('time (hours post signal induction)');
ylabel ('density');
legend ([g1 g2 g3, gtf], {'all data', 'PRG', 'SRG', 'TFs'});

saveas(Fig2b,'6_results\figures\Figure2b.png')


%% Figure3 - predictions of the Uhlitz CYHX dataset
load('4_processed_data\uhl_main\uhl_pred_CYHX.mat');

Fig3 = figure;
set(gcf, 'Position',  [100, 100, 600, 900])

fy = [CYHX_log2fc(:,1) CYHX_log2fc(:,1) ...
      CYHX_log2fc(:,2) CYHX_log2fc(:,2) ...
      CYHX_log2fc(:,3) CYHX_log2fc(:,3)];
fx = [pred_exp(:,1) pred_null(:,1) ...
      pred_exp(:,2) pred_null(:,2) ...
      pred_exp(:,3) pred_null(:,3)];
ftitle = ["CP model - 1h" "null model - 1h" ...
          "CP model - 2h" "null model - 2h" ...
          "CP model - 4h" "null model - 4h"];
%NB. all rsq_null values are less than 0, so take 0
ftext = [round(rsq(1),2) 0 ...
         round(rsq(2),2) 0 ...
         round(rsq(3),2) 0]; 
for i = 1:6
    subplot(3,2,i);
    scatter(fx(:,i), fy(:,i));
    hold on
    plot([-2 8],[-2 8]);
    
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
    xlim([-2 8])
    ylim([-2 8])    
    title(ftitle(i))
    text(3.5, 7, "R2 = " + sprintfc('%0.2f',ftext(i)))
    xlab = text(1.2, -1.5, "predicted expression");
    ylab = text(-1.8, 1, "measured expression");
    set(ylab, 'Rotation', 90)
end

saveas(Fig3,'6_results\figures\Figure3.png')

%% Figure 4a - Qu 2018 CP density plots, GO analysis, and predictions
load('4_processed_data\qu_main\qu_main.mat');
Fig4a = figure;
hold on

[f,xi] = drawDensity(out_opt, qu_days, 0.1);
f = f*0.55;
plot(xi,f, "Color", 'black',"LineWidth",2.0);

daspect([3.5 1 1]) % 0.5*xlim, 1*ylim, and 1
ylim([0 1])
xlim([0 7])
xlabel ('time (days post induction of differentiation)');
ylabel ('density');

saveas(Fig4a,'6_results\figures\Figure4a.png')

%% Figure 4b - EEC syndrome predictions

load('4_processed_data\qu_main\qu_pred_EEC.mat');
Fig4b = figure;
hold on

all_data = [pred_exp, expr_EEC];
all_data = all_data(sum(isnan(all_data),2)==0,:);
[pc, ~, pcvar] = pca (all_data);
pcvar = round(pcvar/sum(pcvar)*100,1);

plot(pc(1:3,1), pc(1:3,2), 'o', 'LineWidth',2, 'MarkerSize',10);
plot(pc(4:6,1), pc(4:6,2), 'o', 'LineWidth',2, 'MarkerSize',10);
text(pc(:,1)+.02, pc(:,2), {'k=1', 'k=2', 'k=3', 'R204W', 'R279H', 'R304W'},...
    'HorizontalAlignment', 'left');

xlabel(['PC 1 (', num2str(pcvar(1)), '%)']);
ylabel(['PC 2 (', num2str(pcvar(2)), '%)']);
grid on;
daspect([0.5 1 1]) % 0.5*xlim, 1*ylim, and 1

saveas(Fig4b,'6_results\figures\Figure4b.png')

end


% local
function P = graph_points(n, cp)
    % generate points from CPs for graphing
    % n = number of CPs
    
    n1 = n-1;

    % calculate the binomial coefficient (x!/(y!(x-y)!)) 
    for i = 0:1:n1
        binomial_coefficient(i+1) = factorial(n1)/(factorial(i)*factorial(n1-i));   
    end

    basis = [];
    UB = [];
    for tau = 0:0.01:1
        for k = 1:n
            UB(k) = binomial_coefficient(k)*((1-tau)^(n-k))*(tau^(k-1));
        end
        basis = cat(1,basis,UB);    
    end

    P = basis*cp;

end

function [f, xi] = drawDensity(out_opt, timeRange, bandwidth)

plotdata = [];

for i = 1:height(out_opt)
    nIS = cell2mat(out_opt(i,2))-2;
    if ~isempty(nIS) && nIS > 0
        geneIS = cell2mat(out_opt(i,3));
            idxnz = find (geneIS (:,2) > 0);
            if ~isempty (idxnz) && geneIS (idxnz (1,1),1) < timeRange (end)
                geneIS = geneIS(idxnz (1,1), 1);
                plotdata = [plotdata; geneIS];
            end
    end
end

%make plot
[f,xi] = ksdensity(plotdata, 'Support', [timeRange(1) timeRange(end)], 'Bandwidth', bandwidth);
end
