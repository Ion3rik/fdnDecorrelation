% Histograms for different mixing configuration
% Authors: Jon Fagerström, Sebastian J. Schlecht
% Updated: 29.5.2023

clear; clc;
%% INIT
rng(5)
numInstances = 100; % average median estimate over numInstances
fs = 48000;

mixingType = 'orthogonal';

% define FDN
N = 4;

maxIODelay = 100;

%% CORRELATION ANALYSIS

for instance = 1:numInstances
    A = fdnMatrixGallery(N,'orthogonal');       % feedback matrix
    m = randi([300,10000],[1,N]);                % FDN delays
    MixIn = fdnMatrixGallery(N,mixingType);   % input mixing matrix
    MixOut = fdnMatrixGallery(N,mixingType);  % output mixing matrix
    % Construct P matrix
    P = loopTF(m,A);
    adjMat = adjPoly(P,'z^1');

    for delayConfigs = 1:4
        inputDelays = randi([1,maxIODelay],[1,N]);
        outputDelays = randi([1,maxIODelay],[1,N]);

        switch delayConfigs
            case 1 % no delays
                Btemp = eye(N);
                Ctemp = eye(N);
            case 2 % input delays
                Btemp = constructDelayMatrix(inputDelays);
                Ctemp = eye(N);
            case 3 % output delays
                Btemp = eye(N);
                Ctemp = constructDelayMatrix(outputDelays);
            case 4 % both delays
                Btemp = constructDelayMatrix(inputDelays);
                Ctemp = constructDelayMatrix(outputDelays);
        end
        for mixingConfigs = 1:4
            switch mixingConfigs
                case 1 % no  mixing
                    B = Btemp;
                    C = Ctemp;
                case 2 % input mixing
                    B = matrixConvolution(Btemp, MixIn);
                    C = Ctemp;
                case 3 % output mixing
                    B = Btemp;
                    C = matrixConvolution(MixOut, Ctemp);
                case 4 % both mixing
                    B = matrixConvolution(Btemp, MixIn);
                    C = matrixConvolution(MixOut, Ctemp);
            end

            feedforwardPath = matrixConvolution(C, matrixConvolution(adjMat,B));
            %feedforwardPath = adjMat;
            %% Compute correlation
            maxLag = sum(m);

            % adjMat
            %[~, ~, ~, maxCorrelation] = corr(feedforwardPath, maxLag, [], []);
            x = maxCorr(feedforwardPath);
            % filter upper triangle
            x = triu(x,1);
            xx = x(:);
            xx( abs(xx) < eps ) = [];
            maxCorrelation{mixingConfigs,delayConfigs,instance} = abs(xx);
        end
    end
end

%% PLOT
fig7 = figure('units','centimeters','position',[0,0,30,10]);
hold on;
t = tiledlayout(4,4);
columntitle = {'No delays','Input delays','Output delays','Both delayed'};
rowtitle = {'No mixing','Input mixing','Output mixing','Both mixed'};
for mixingConfigs = 1:4
    for delayConfigs = 1:4
        allVal = [];
        for instance = 1:numInstances
            allVal = [allVal ,maxCorrelation{mixingConfigs,delayConfigs,instance}];
        end
        edges = linspace(0,1,40);
        nexttile;
        histogram(allVal(:),edges,'Normalization','probability');
        ylim([0 0.1])
    end
end
% Specify common title, X and Y labels
% title(t, 'Common title')
xlabel(t, 'Absolute Correlation')
ylabel(t, 'Probability')

for mixingConfigs = 1:4
     annotation('textbox', [0,1.15-mixingConfigs*0.22, 0, 0], 'string', rowtitle(mixingConfigs))
end

for delayConfigs = 1:4
     annotation('textbox', [-0.05+delayConfigs*0.21, 0.99, 0.1, 0.], 'string', columntitle(delayConfigs),'EdgeColor','none')
end
