% startup
clear; clc; close all;

addpath(genpath('./auxiliary'))
addpath(genpath('./Externals/fdnToolbox'))

%% List File and Product dependencies
% files = dir('./**/*.m');
% ff = {files.name}'
% [fList,pList] = matlab.codetools.requiredFilesAndProducts(ff)
% pList.Name