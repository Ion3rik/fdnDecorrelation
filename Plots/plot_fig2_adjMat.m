% Plot the Adjugate %
% Author: Jon Fagerström %
% Updated: 29.5.2023 %


% Figure 2: Time-domain filter coefficients of the adjugate adj(P(z)) of a 
% MIMO FDN with four delay lines ($\matSize = 4$) and \edit{with} a random 
% orthogonal feedback matrix A. The delays are m = [977, 683, 981, 801] samples. 
% The entire adjugate adj(P(z)) is displayed without truncation. 
% Only the non-zero values are drawn with stems for better readability.

clear; clc;
%% INIT
rng(4) % init random number generator
fs = 48000; % sampling rate
N = 4; % FDN size

%% COMPUTE THE ADJACENT MATRIX

numInput = N;
numOutput = N;
B = eye(N,numInput);
C = eye(numOutput,N);

D = zeros(numOutput,numInput);
m = randi([300,1000],[1,N]);
A = fdnMatrixGallery(N,'orthogonal');

% Construct P matrix
P = loopTF(m,A);
adjMat = adjPoly(P,'z^1');

%% PLOT
% plot only non-zeros
eps = 10^(-6);
adjMat(abs(adjMat)<eps) = nan;

c1 = [0 0.4470 0.7410];
fig2 = figure;
h = plotImpulseResponseMatrix([], adjMat, 'LineWidth', 1.5, 'MarkerSize', 3);
han = axes(fig2,'visible','off'); han.XLabel.Visible = 'on'; han.YLabel.Visible = 'on';
set(h, 'xlim',[-10, size(adjMat,3)+100], 'ylim', [-1 1]);
xlabel(han,'Time (samples)'); ylabel(han,'Sample value');
fig2.Position = [100 100 900 350];
ax = findobj(gcf,'type','axes');
for i = 3:18
    yline(ax(i),0, 'Linewidth', 1, 'color', c1)
end

