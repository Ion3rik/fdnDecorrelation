% Generate input/output delay filter results
% TABLE III: Median correlation metric Φ for different delay filter configurations
% Jon Fagerström
% 8.6.2022

clear; clc; close all;
rng(5)
numInstances = 10; % average median estimate over numInstances
fs = 48000;

simo = false; % toggle SIMO mode
mixingType = 'orthogonal';

% define FDN
NN = [4 8 16];


maxIODelay = 100;

for instance = 1:numInstances
    instance
    for sizes = 1:numel(NN)
        N = NN(sizes);
        A = fdnMatrixGallery(N,'orthogonal');       % feedback matrix
        MixIn = fdnMatrixGallery(N,mixingType);   % input mixing matrix
        MixOut = fdnMatrixGallery(N,mixingType);  % output mixing matrix
        m = randi([300,10000],[1,N]);
        inputDelays = randi([1,maxIODelay],[1,N]);
        outputDelays = randi([1,maxIODelay],[1,N]);
        for cases = 1:4
            cases
            switch cases
                case 1 % vanilla
                    B = eye(N);
                    C = eye(N);
                case 2 % input
                    B = constructDelayMatrix(inputDelays);
                    B = matrixConvolution(B, MixIn);
                    C = eye(N);
                case 3 % output
                    B = eye(N);
                    C = constructDelayMatrix(outputDelays);
                    C = matrixConvolution(MixOut, C);
                case 4 % both
                    B = constructDelayMatrix(inputDelays);
                    B = matrixConvolution(B, MixIn);
                    C = constructDelayMatrix(outputDelays);
                    C = matrixConvolution(MixOut, C);
            end
                % Construct P matrix
                P = loopTF(m,A);
                adjMat = adjPoly(P,'z^1');
                feedforwardPath = matrixConvolution(C, matrixConvolution(adjMat,B));

                % Compute correlation
                maxLag = sum(m);
                %[~, ~, ~, maxCorrelation] = corr(feedforwardPath, maxLag, [], []);
                maxCorrelation = maxCorr(feedforwardPath);

                % filter upper triangle
                x = triu(maxCorrelation,1);
                xx = x(:);
                xx( abs(xx) < eps ) = [];
                medCurrent{cases,sizes,instance} = abs(xx);
        end
    end
end

%% Collect data
for sizes = 1:numel(NN)
    for cases = 1:4
        allVal = [];
        for instance = 1:numInstances
            allVal = [allVal medCurrent{cases,sizes,instance}];
        end
        med(cases,sizes) = median(allVal(:));
        IQR(cases,sizes) = iqr(allVal(:));
    end
end


filename = ['table3.mat'];
save(filename);

% printMatLatex(med,'format','%1.4f')  % print
% printMatLatex(IQR,'format','%1.4f')