% Delays + Mixing median correlation
% Jon Fagerström
% TABLE II: Median correlation metric Φ for different mixing configurations.
% 24.5.2022

clear; clc; close all;
rng(5)
numInstances = 10; % average median estimate over numInstances
fs = 48000;
impulseResponseLength = fs/10;

simo = false; % toggle SIMO mode
mixingType = 'orthogonal';

% define FDN
N = 4;

maxIODelay = 100;

for instance = 1:numInstances
    A = fdnMatrixGallery(N,'orthogonal');       % feedback matrix
    %A = constructVelvetFeedbackMatrix(N,3,3);
    m = randi([300,10000],[1,N]);                % FDN delays
    MixIn = fdnMatrixGallery(N,mixingType);   % input mixing matrix
    MixOut = fdnMatrixGallery(N,mixingType);  % output mixing matrix
    % Construct P matrix
    P = loopTF(m,A);
    adjMat = adjPoly(P,'z^1');

    for delayConfigs = 1:4
        delayConfigs
        %         Btemp = zeros(N, N, maxIODelay);
        %         Ctemp = zeros(N, N, maxIODelay);
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
            maxCorrelation = maxCorr(feedforwardPath);
            % filter upper triangle
            x = triu(maxCorrelation,1);
            xx = x(:);
            xx( abs(xx) < eps ) = [];
            medCurrent{mixingConfigs,delayConfigs,instance} = abs(xx);
        end
    end
end

%% Collect data
for delayConfigs = 1:4
    for mixingConfigs = 1:4
        allVal = [];
        for instance = 1:numInstances
            allVal = [allVal medCurrent{mixingConfigs,delayConfigs,instance}];
        end
        med(mixingConfigs,delayConfigs) = median(allVal(:));
        IQR(mixingConfigs,delayConfigs) = iqr(allVal(:));
    end
end


filename = ['table2.mat'];
save(filename);

% printMatLatex(med,'format','%1.4f')
% printMatLatex(IQR,'format','%1.4f')