% Plot MIMO and SIMO adjacent matrix
% Authors: Sebastian J. Schlecht, Jon Fagerström
% Updated: 29.5.2023

clear; clc;

%% INIT
rng(14)
fs = 48000;
mixingType = 'orthogonal';

%% FDN PARAMS
N = 4;

maxIODelay = 100;

% scattering matrix params
sparsity = 17;
numStages = 2;
A = constructVelvetFeedbackMatrix(N,numStages,sparsity);

m = randi([300,1000],[1,N]);                % FDN delays
P = loopTF(m,A);
adjMat = adjPoly(P,'z^1');
adjMat = adjMat(:,:,1:2000-1);

inputDelays = randi([1,maxIODelay],[1,N]);
B = constructDelayMatrix(inputDelays);
B = sum(B,2);

%% PLOT

c1 = [0 0.4470 0.7410];
fig8a = figure('Units','pixels');
fig8a.Position = [100 100 800 450];

h = plotImpulseResponseMatrix([], removeZeros(adjMat), 'LineWidth', 1.5, 'MarkerSize', 3);
han = axes(fig8a,'visible','off'); han.XLabel.Visible = 'on'; han.YLabel.Visible = 'on';
set(h, 'xlim',[-10, size(adjMat,3)+100], 'ylim', [-1 1], 'FontSize', 14);
xlabel(han,'Time (samples)'); ylabel(han,'Sample value');
ax = gca;
ax.FontSize = 14;

ax = findobj(gcf,'type','axes');
for i = 3:18
    yline(ax(i),0, 'Linewidth', 1, 'color', c1)
end

fig8b = figure('Units','pixels');
fig8b.Position = [900 100 280 450];
h = plotImpulseResponseMatrix([], removeZeros(sum(adjMat,2)), 'LineWidth', 1.5, 'MarkerSize', 3);
han = axes(fig8b,'visible','off'); han.XLabel.Visible = 'on'; han.YLabel.Visible = 'on';
set(h, 'xlim',[-10, size(adjMat,3)+100], 'ylim', [-1 1], 'FontSize', 14);
xlabel(han,'Time (samples)'); ylabel(han,'Sample value');
ax = gca;
ax.FontSize = 14;
ax = findobj(gcf,'type','axes');
for i = 3:6
    yline(ax(i),0, 'Linewidth', 1, 'color', c1)
end


%% FUNCTIONS
function tSettings(fig, font)
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Sample value');
    xlabel(han,'Time (samples)');
    ax = gca;
    ax.FontSize = font;
end

function A = removeZeros(A)
% remove zeros
eps = 10^(-6);
A(abs(A)<eps) = nan;
end