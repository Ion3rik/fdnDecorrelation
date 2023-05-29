% Histograms of the Inter-channel Max Correlation for different Feedback
% matrices
% Authors: Jon Fagerström, Sebastian J. Schelcht
% Updated: 29.5.2023

clear; clc;
%% INIT
rng(5)
numInstances = 10; % number of randomized parameter inits (fb-matrix and delay lengths)
fs = 48000; % sampling rate

% define FDN
N = 4;

% scattering matrix params
sparsity = 3;
numStages = 3;

%% CORRELATION ANALYSIS
numCases = 6;
for instance = 1:numInstances
    m = randi([300,10000],[1,N]);
    for cases = 1:numCases
        switch cases
            case 1
                A = fdnMatrixGallery(N,'orthogonal');
            case 2
                A = fdnMatrixGallery(N,'Hadamard');
            case 3
                A = fdnMatrixGallery(N,'Householder');
            case 4
                A = fdnMatrixGallery(N,'circulant');
            case 5
                A = constructVelvetFeedbackMatrix(N,numStages,sparsity);
            case 6
                A = constructCascadedParaunitaryMatrix(N,numStages);
        end

        % Construct P matrix
        P = loopTF(m,A);
        adjMat = adjPoly(P,'z^1');
        
        %% Compute maximum correlation
        x = maxCorr(adjMat);
        % filter upper triangle
        x = triu(x,1);
        xx = x(:);
        xx( abs(xx) < eps ) = [];
        maxCorrelation{cases,instance} = abs(xx);
    end
end
%% PLOT 
close all;
set(0,'defaultfigurecolor',[1 1 1])
lineStyle{1} = '-'; lineStyle{2} = '--'; lineStyle{3} = ':';
lineStyle{4} = '-'; lineStyle{5} = '--'; lineStyle{6} = ':';
fig5 = figure('units','centimeters','position',[0,0,30,10]);
t = tiledlayout('flow'); 
nexttile; hold on;
xlabel(t, 'Absolute Correlation')
ylabel(t, 'Probability')
for cases = 1:numCases
    allVal = [];
    for instance = 1:numInstances
        allVal = [allVal maxCorrelation{cases,instance}]; % collect all values
    end

    if(cases == 4)
        legend('Random Orthogonal','Hadamard','Householder')
        hold off;
        nexttile; hold on;
        set(gca,'ColorOrderIndex',4)
        
    end

    edges = linspace(0,1,(N*3)^2);
    histogram(allVal(:),edges,'Normalization','probability', 'LineStyle',lineStyle{cases});
    ylim([0 0.2]);
end
legend('Circulant','Velvet Scattering','Dense Scattering')
