% PSHA calculation considering two ruptures.
% The calculation follows the example of Section 6.3.2
%
% Created by Jack Baker

clear; close all; clc
addpath('../utils/')

% basic setup
x = logspace(log10(0.001), log10(2), 100); % IM values to consider 
T = 1; % period of interest
IM_label = 'SA(1 s)';
gmpeFlag = 1; % =1 for BJF97, =2 for CY14

% specify colors and line styles for plots
colorspec{1} = [56 95 150]/255;
colorspec{2} = [207 89 33]/255;
colorspec{3} = [158 184 219]/255;


% seismicity parameters
Fault_Type = 1; % 1 is strike slip
Vs30 = 500;

% plotting parameters
figureAxisLimits = [0.05 max(x) 1e-5 1e-1];
figureXTickVals = [0.05 0.1 0.5 1 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first rupture 

lambda_A = 1/100;
M_A = 6.5;
R_A = 10;

% compute rates (and intermediate results) for specific IM levels
[medianIM, sigmaIM]  = gmpe_bjf97(M_A, R_A, T, Fault_Type, Vs30);
imLevel(1) = 0.2;
imLevel(2) = 0.5;
imProbabilitiesA = 1 - normcdf(log(imLevel),log(medianIM),sigmaIM)
imRateA = lambda_A * imProbabilitiesA; % get rates for two example cases

% compute rates for a range of IM levels
p_A = 1 - normcdf(log(x),log(medianIM),sigmaIM);
lambda_IM_A = lambda_A * p_A; % IM rates from rup_1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% second rupture

% define second fault
lambda_B = 1/500;
M_B = 7.5;
R_B = 10;

% compute rates (and intermediate results) for specific IM levels
[medianIM, sigmaIM]  = gmpe_bjf97(M_B, R_B, T, Fault_Type, Vs30);
imProbabilitiesB = 1 - normcdf(log(imLevel),log(medianIM),sigmaIM)
imRateB = lambda_B * imProbabilitiesB; % get rates for two example cases
imRateTot = imRateA + imRateB; 

% compute rates for a range of IM levels
p_B = 1 - normcdf(log(x),log(medianIM),sigmaIM);
lambda_IM_B = lambda_B * p_B; % IM rates from rup_2

lambda_IM_Tot = lambda_IM_A + lambda_IM_B;

figure
loglog(x, lambda_IM_Tot,'-', 'linewidth', 2, 'color', colorspec{1})
hold on
loglog(x, lambda_IM_A,'-',   'linewidth', 2, 'Color', colorspec{2})
loglog(x, lambda_IM_B,'-',   'linewidth', 2, 'Color', colorspec{3})

plot(imLevel, imRateTot, 'o', 'MarkerEdgeColor', colorspec{1}) 
plot(imLevel, imRateA,   'o', 'MarkerEdgeColor', colorspec{2}) 
plot(imLevel, imRateB,   'o', 'MarkerEdgeColor', colorspec{3}) 

% annotate text results for example cases
text1 = ['\lambda(' IM_label ' > ' num2str(imLevel(1)) ' g) = ' num2str(imRateTot(1),3)];
text2 = ['\lambda(' IM_label ' > ' num2str(imLevel(2)) ' g) = ' num2str(imRateTot(2),3)];
text(imLevel(1)*1.1, imRateTot(1)*1.1,text1)
text(imLevel(2)*1.05, imRateTot(2)*1.2,text2)

xlabel(['Spectral Acceleration, ' IM_label ' [g]'])
ylabel('Annual rate of exceedance, \lambda')
axis(figureAxisLimits)
legend('Total hazard', 'rup_1', 'rup_2');





