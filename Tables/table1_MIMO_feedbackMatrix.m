% Example on the decorrelation depending on feedback matrix
% TABLE 1:  Median correlation metric Φ for different feedback matrix types and sizes.
% 7.5.2021
% Jon Fagerström
% Sebastian J. Schlecht, Tuesday, 17. May 2022
clear; close all;
rng(5)
numInstances = 10; % average median estimate over numInstances
fs = 48000;
impulseResponseLength = fs/10;

% define FDN
NN = [4 8 16];
%NN = 4;

% scattering matrix params
sparsity = 3;
numStages = 3;

numCases = 6;
for instance = 1:numInstances
    instance
    for sizes = 1:numel(NN)
        N = NN(sizes);
        m = randi([300,10000],[1,N]);
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
            
            %% Compute correlation
            maxLag = sum(m);

            maxCorrelation = maxCorr(adjMat);
            % filter upper triangle
            x = triu(maxCorrelation,1);
            xx = x(:);
            xx( abs(xx) < eps ) = [];
            medCurrent{cases,sizes,instance} = abs(xx);
        end
    end
end

%% Collect data
filename = ['table1.mat'];
% load(filename)
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
save(filename);

% printMatLatex(med,'format','%1.4f')
% printMatLatex(IQR,'format','%1.4f')