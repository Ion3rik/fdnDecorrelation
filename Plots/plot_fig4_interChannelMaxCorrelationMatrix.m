% Inter-channel Max Correlation
% Authors: Jon Fagerström, Sebastian J. Schelcht
% Updated: 29.5.2023

clear; clc;
%% INIT
rng(4) % init random number generator
fs = 48000; % sampling rate
impulseResponseLength = fs/10;
N = 4; % FDN Size

%% PARAMTERIZE FDN
numInput = N;
numOutput = N;
B = eye(N,numInput);
C = eye(numOutput,N);

D = zeros(numOutput,numInput);
m = randi([300,1000],[1,N]);
A = fdnMatrixGallery(N,'orthogonal');
        
%% DECOMPOSE
% Construct P matrix
P = loopTF(m,A);
adjMat = adjPoly(P,'z^1');

%% COMPUTE INTER-CHANNEL MAX CORRELATION
maxCorrelation = maxCorr(adjMat);
%% PLOT
fig4 = figure;
fig4.Position = [0 0 700 600];
plotHeatMap(abs(maxCorrelation), [0 1])
xlabel('ij'); ylabel('kl')


%% FUNCTIONS
function h = plotHeatMap(matrix, limits)
    dim = size(matrix,1);
    N = sqrt(dim);
    for i = 1:dim
        [a,b] = ind2sub(N, i);
        coordLabel{i} = [num2str(a) num2str(b)];
    end
    h = heatmap(coordLabel, coordLabel,matrix);
    grid on;
    colormap(gray)
    h.ColorLimits = limits; 
    Ax = gca;
end
