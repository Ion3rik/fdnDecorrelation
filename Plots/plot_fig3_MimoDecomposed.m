% MIMO FDN Decomposed Impulse Response
% Authors: Sebastian J. Schlecht, Jon Fagerström %
% Updated: 29.5.2023 %

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
p = generalCharPoly(m,A);

%% COMPUTE IMPULSE RESPONSES
impulseResponse = permute( dss2impz(impulseResponseLength, m, A, eye(N), eye(N), D), [2 3 1]);
impulseResponse = matrixConvolution(C, matrixConvolution(impulseResponse,B));

undrivenResponse = impz(1,p,impulseResponseLength);
undrivenResponse = [zeros(min(m), 1); undrivenResponse];
driveResponse =  matrixConvolution(C, matrixConvolution(adjMat,B));
driveUndrive = conv(undrivenResponse,squeeze(driveResponse(1,1,:)));

% check that the results are the same
max(abs(squeeze(impulseResponse(1,1,:)) - driveUndrive(1:length(impulseResponse))))



%% PLOT
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.9290 0.6940 0.1250];
mSize = 3;

% Remove zeros for plotting
eps = 10^(-6);
fir = squeeze(driveResponse(1,1,:)); fir(abs(fir)<eps) = nan;
iir = undrivenResponse; iir(abs(iir)<eps) = nan;
ir = driveUndrive; ir(abs(ir)<eps) = nan;


fig3 = figure; 
x = axes;
grid on;
ylim([-1 5]);

hold on

a = stem(fir+4, 'filled',  'Linewidth', 2, 'color', c1, 'BaseValue',4, 'MarkerSize', mSize);
al = a.BaseLine; al.Color = 'none';
xlabel('Time (samples)')
ylabel('Sample value')
plot(1:length(fir),4*ones(size(fir)), 'color', c1, 'Linewidth', 1);
xlim([1 impulseResponseLength]); ylim([-1 5]); 
y = axes('position', x.Position);
set(y, 'Color', 'none','xtick',[],'ytick',[]); 
hold on
b = stem(iir+2, 'filled',  'Linewidth', 2, 'color', c2, 'BaseValue',2, 'MarkerSize', mSize);
bl = b.BaseLine; bl.Color = 'none';
plot(1:length(iir),2*ones(size(iir)), 'color', c2, 'Linewidth', 1);
xlim([1 impulseResponseLength]); ylim([-1 5]); 
z = axes('position', x.Position);
set(z, 'Color', 'none','xtick',[],'ytick',[]); 
xlim([1 impulseResponseLength]); ylim([-1 5]);
hold on
c = stem(ir, 'filled',  'Linewidth', 2, 'color', c3, 'BaseValue',0, 'MarkerSize', mSize);
cl = c.BaseLine; cl.Color = 'none';
plot(1:length(ir),zeros(size(ir)), 'color', c3, 'Linewidth', 1);
legend([a; b; c],{'Feedforward','Recursive','Impulse response'})


fig3.Position = [100 100 900 280];
