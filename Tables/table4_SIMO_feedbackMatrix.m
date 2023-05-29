% SIMO Median correlation
% TABLE IV: Median correlation metric Φ for various SIMO FDNs
% Jon Fagerström
% 14.6.2022

clear; clc; close all;
rng(7)

numInstances = 10; % average median estimate over numInstances
fs = 48000;

mixingType = 'orthogonal';

% define FDN
NN = [4 8 16];
% NN = 4;

% scattering matrix params
sparsity = 8;
numStages = 3;

maxIODelay = 100;

numCases = 6;

for instance = 1:numInstances
    instance
    for sizes = 1:numel(NN)
        N = NN(sizes);
        MixIn = fdnMatrixGallery(N,mixingType);   % input mixing matrix
        MixOut = fdnMatrixGallery(N,mixingType);  % output mixing matrix
        m = randi([300,10000],[1,N]);
        inputDelays = randi([1,maxIODelay],[1,N]);
        outputDelays = randi([1,maxIODelay],[1,N]);

        B = ones(N,1);
        C = eye(N);

%         B = constructDelayMatrix(inputDelays);
%         B = matrixConvolution(B, MixIn);
%         B = sum(B,2);
%         C = constructDelayMatrix(outputDelays);
%         C = matrixConvolution(MixOut, C);

        for cases = 1:numCases
            cases
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

                
%                 C = matrixConvolution(A,constructDelayMatrix(outputDelays));
% 
%                 feedforwardPath = matrixConvolution(C, matrixConvolution(adjMat,B));
                feedforwardPath = matrixConvolution(adjMat,B);
                        
                %% Compute correlation
                maxLag = sum(m);

                % adjMat
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
    for cases = 1:numCases
        allVal = [];
        for instance = 1:numInstances
            allVal = [allVal medCurrent{cases,sizes,instance}];
        end
        med(cases,sizes) = median(allVal(:));
        IQR(cases,sizes) = iqr(allVal(:));
    end
end

filename = ['table4.mat'];
save(filename);

% printMatLatex(med,'format','%1.4f')
% printMatLatex(IQR,'format','%1.4f')