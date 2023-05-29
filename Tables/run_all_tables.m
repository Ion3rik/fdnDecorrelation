%% Run all table scripts
tic
table1_MIMO_feedbackMatrix
table2_mixing
table3_delayFilter
table4_SIMO_feedbackMatrix
toc

%% Print all tables

load('table1.mat');
printTable({'Random Orthogonal','Hadamard','Householder','Circulant','Velvet Scattering','Dense Scattering'}, combineMedAndIQR(med,IQR),'%1.3f');

load('table2.mat');
printTable({'No mixing','Input mixing','Output mixing','Both mixed'}, combineMedAndIQR(med,IQR),'%1.3f');

load('table3.mat');
printTable({'No filters','Input filters','Output filters','Both filters'}, combineMedAndIQR(med,IQR),'%1.3f');

load('table4.mat');
printTable({'Random Orthogonal','Hadamard','Householder','Circulant','Velvet Scattering','Dense Scattering'}, combineMedAndIQR(med,IQR),'%1.3f');
 

%% FUNCTIONS

function printTable(names, mat,format)

[m,n] = size(mat);
ff = [' %s & ', repmat([format ' & '],1,n-1), [' ' format ' \\\\ \n']];
C = [names.' num2cell(mat)].';
fprintf(ff,C{:})
fprintf('\n')
end

function out = combineMedAndIQR(med, IQR)
    [m,n] = size(med);
    out = zeros(m,2*n);
    out(:,1:2:end) = med;
    out(:,2:2:end) = IQR;
end