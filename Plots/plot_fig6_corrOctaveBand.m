% Compute SIMO FDN median correlation metric at octave bands
% Authors: Jon Fagerström 
% Updated: 29.5.2023

clear; clc;
%% INIT
rng(5)
numInstances = 10; % average median estimate over numInstances
fs = 48000; % sampling rate

%% PARAMETERS
NN = 4; % FDN size

% scattering matrix params
sparsity = 3;
numStages = 3;

% Octave Filter Params
fc = [63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000];
numBands = numel(fc);

for band = 1:numBands
    octFilt{band} = octaveFilter(fc(band),"1 octave");
end

for instance = 1:numInstances
    for sizes = 1:numel(NN)
        N = NN(sizes);
        m = randi([300,1000],[1,N]);
        A = constructVelvetFeedbackMatrix(N,numStages,sparsity);

        % Construct P matrix
        P = loopTF(m,A);
        adjMat = adjPoly(P,'z^1');
        % Octave Bands
        for band = 1:numBands+1
            if band > numBands
                adjMatBand = adjMat; % broadband
            else
                adjMatBand = octaveFilterMat(adjMat,octFilt{band}); % filter
            end
            % Compute correlation
            maxLag = sum(m);
            x= maxCorr(adjMatBand);
            % filter upper triangle
            x = triu(x,1);
            xx = x(:);
            xx( abs(xx) < eps ) = [];
            maxCorrelation{band,sizes,instance} = abs(xx);
        end
    end
end

%% Collect data

for sizes = 1:numel(NN)
    for bands = 1:numBands+1
        allVal = [];
        for instance = 1:numInstances
            allVal = [allVal maxCorrelation{bands,sizes,instance}];
        end
        med(bands,sizes) = median(allVal(:));
    end
end

%% PLOT
h = 200; w = 400;
fig6 = figure('Renderer', 'painters', 'Position', [1000 310 w h]);
for sizes = 1:numel(NN)
    plot(fc, med(1:end-1,sizes), 'LineWidth', 2); hold on;
end
legend(string(NN));
xlabel('Frequency [Hz]')

ylabel('Absolute Correlation')
xlim([fc(1), fc(end)])
ax = gca;
ax.FontSize = 12;
set(gca,'XScale', 'log', 'XTick',[50 100 250 500 1000 2000 4000 8000 16000], 'XTicklabel',{'50', '100', '250', '500', '1k', '2k', '4k', '8k', '16k'});
grid on;

%% FUNCTIONS

function irFilt = octaveFilterMat(ir, octFilt)
    [n,m,~] = size(ir); % dimensions
    irFilt = zeros(size(ir));
    for i = 1:n
        for j = 1:m
            irFilt(i,j,:) = octFilt(squeeze(ir(i,j,:)));
        end
    end
end
