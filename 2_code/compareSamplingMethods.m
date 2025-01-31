% two demo cases comparing our weighted bootstrap method, spline
% interpolation, and Gaussian process (GP) regression.
% Rosemary Yu. Last update: 2025-01-31

clear 
clc
close all

fig = figure;
tiledlayout(1,2,'TileSpacing','compact')
nexttile

x = 0:5; 
y = [1 1 1 -1 -1 -1]; 
drawComparisons(x,y)
nexttile

x = 0:5; 
y = [0 0 1 0 0 0]; 
drawComparisons(x,y)


% local
function drawComparisons(x, y)
    % compare bootstrapping, spline, and GP regression 
    
    % bootstrap
    tau = getTau(horzcat(x.', y.'));
    b = bootTimeSeries(x,y, 100, tau);

    % spline
    xpred = b(:,1);
    ypred_spline = spline(x,y,xpred.');

    % GP regression with adjusted parameters 
    % (with default parameters, the model breaks in the 2nd demo case)
    gpr2 = fitrgp(x.', y.', 'FitMethod','none', 'Sigma', 0.2);
    [ypred_gpr,~,~] = predict(gpr2,xpred);

    %plot all
    plot(x,y,'o',...                  % data
        xpred, b(:,2), '--',...       % bootstrap
        xpred,ypred_spline,'-', ...   % spline
        xpred,ypred_gpr,'-.')         % GP regression
    
    legend('data','bootstrap','spline','GP regression')
end
