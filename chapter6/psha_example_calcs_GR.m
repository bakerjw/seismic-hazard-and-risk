% PSHA calculation using a Gutenberg-Richter magnitude distribution.
% The calculation follows the example of Section 6.3.3
%
% Created by Jack Baker


clear; close all; clc
addpath('../utils/')



% basic setup
x = logspace(log10(0.001), log10(2), 100); % IM values to consider 
T = 1; % period of interest
IM_label = 'SA(1 s)';
Axis_label = 'Spectral Acceleration, SA(1 s) [g]';
gmpeFlag = 1; % =1 for BJF97, =2 for CY14

% specify colors and line styles for plots
colorspec{1} = [56 95 150]/255;
colorspec{2} = [207 89 33]/255;
colorspec{3} = [158 184 219]/255;


% seismicity parameters
R = 10;
Rrup = R; Rjb = R; 
rup.Fault_Type = 1; % 1 is strike slip
rup.Vs30 = 500;
rup.R = 10;
% CY parameters
rup.Ztor  = 0;
rup.delta = 90;
rup.rupLambda = 0;
rup.Z10 = 999;
rup.Fhw = 0;
rup.FVS30 = 0;
rup.region = 1;


% plotting parameters
figureAxisLimits = [0.05 max(x) 1\0.99e-5 1e-1];
figureXTickVals = [0.05 0.1 0.5 1 2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% point source with G-R magnitudes

M = [5 6 7 8];     

% From Table 3.5, fixed rate of M>5, M_max = 8
lambda_M = [0.05	0.03153	0.01988	0.01252	0.007882	0.004955	0.003108	0.001942	0.001207	0.0007432	0.0004505	0.0002657	0.0001492	7.57E-05	2.93E-05];
M_vals = [5.1	5.3	5.5	5.7	5.9	6.1	6.3	6.5	6.7	6.9	7.1	7.3	7.5	7.7	7.9];

x_example = 0.2; % example values for table
[lambda, example_output, disagg] = fn_PSHA_given_M_lambda(lambda_M, M_vals, T, x, x_example, rup, gmpeFlag);

x_example2 = 0.5; % output results for a second threshold
[lambda2, example_output2, disagg2] = fn_PSHA_given_M_lambda(lambda_M, M_vals, T, x, x_example2, rup, gmpeFlag);


%% hazard curve
figure
loglog(x, lambda.x,'-', 'linewidth', 2, 'color', colorspec{1})
hold on
plot(x_example, lambda.example, 'o', 'color', colorspec{1})
plot(x_example2, lambda2.example, 'o', 'color', colorspec{1})

% annotate text results for example cases
text1 = ['\lambda(' IM_label ' > ' num2str(x_example) ' g) = ' num2str(lambda.example,3)];
text2 = {['\lambda(' IM_label ' > ' num2str(x_example2) ' g) ']; ['      = ' num2str(lambda2.example,3)]};

text(x_example*1.1, lambda.example*1.2,text1)
text(x_example2*1.1, lambda2.example*1.2,text2)

xlabel(Axis_label)
ylabel('Annual rate of exceedance, \lambda')
axis(figureAxisLimits)
set(gca, 'xtick', figureXTickVals)

%% output a subset of the hazard curve for use in a table
imSmall = [1e-3 0.1:0.1:1];
ratesSmall = exp(interp1(log(x),log(lambda.x),log(imSmall))); %loglog interpolate
hazTable = [imSmall' ratesSmall'];



%% disaggregation
figure
subplot(1,2,1)
bar(M_vals, disagg.example, 1,'FaceColor', colorspec{3})
hold on
plot( disagg.Mbar*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' > ' num2str(x_example) ' g)']);
axis([5 8 0 0.2])

subplot(1,2,2)
bar(M_vals, disagg2.example, 1,'FaceColor', colorspec{3})
hold on
plot( disagg2.Mbar*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' > ' num2str(x_example2) ' g)']);
axis([5 8 0 0.2])

mBar = [disagg.Mbar disagg2.Mbar] 

% tabulate output
disagg_table = [M_vals' disagg.example' disagg2.example'];



%% Metrics to evaluate calculations and figure

% im with given rate
rateTarg = 1/1000;
imTarg = interp1(ratesSmall, imSmall, rateTarg) % linear interpolation
imTarg = exp(interp1(log(ratesSmall), log(imSmall), log(rateTarg))) % log interpolation

lnImManual  = ( (log (0.2) - log (0.3)) * (log (1E-3) - log (6.81E-4)) ) / (log (2.7E-3) - log (6.81E-4)) + log(0.3) % manual log interpolation
imManual = exp(lnImManual)

% hazard curves slope
imSlope = [0.2 0.3];
rateSlope = exp(interp1(log(x),log(lambda.x),log(imSlope))); %loglog interpolate
kEst = - (log(rateSlope(1)) - log(rateSlope(2)))/ (log(imSlope(1))- log(imSlope(2)));
k0Est = rateSlope(1) / exp(-kEst * log(imSlope(1)));
lambdaPowerLaw = k0Est * exp(-kEst*log(x));

%% hazard curve derivative

dLambda = -diff([ratesSmall 0]);

figure
subplot(1,2,1);
bar(imSmall+0.05, dLambda, 1,'FaceColor', colorspec{3})
axis([0 1 0 0.05])
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel(Axis_label)
ylabel('\Delta \lambda_i ')
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')

% finer discretization
xFine = 0.01:0.01:1;
lambdaFine = exp(interp1(log(x),log(lambda.x),log(xFine))); %loglog interpolate
dLambdaFine = -diff([lambdaFine 0]);

subplot(1,2,2);
bar(xFine+0.005, dLambdaFine, 1,'FaceColor', colorspec{3})
axis([0 1 0 0.008])
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')

xlabel(Axis_label)
ylabel('\Delta \lambda_i ')
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')



%% summary plot
figure
h1 = loglog(x, lambda.x,'-', 'linewidth', 2, 'color', colorspec{1});
hold on
plot(imTarg, rateTarg, 'ok')
h2 = plot(x, lambdaPowerLaw,'-', 'linewidth', 2, 'color', colorspec{3});
plot([0.01 imTarg imTarg], [rateTarg rateTarg 1e-10], ':k', 'linewidth', 1)
% annotate text results for example cases
text1 = ['\lambda(' IM_label ' > ' num2str(x_example) ' g) = ' num2str(lambda.example,3)];
text2 = ['\lambda(' IM_label ' > ' num2str(imTarg,3) ' g) = ' num2str(rateTarg,3)];
% text3 = ['\lambda(M \geq M_{min}) = ' num2str(lambda_M(1))];
text3 = ['\Sigma_i \lambda(rup_i) = ' num2str(lambda_M(1))];
% text(x_example*1.05, lambda.example*1.2,text1)
text(imTarg*1.05, rateTarg*1.5,text2)
text(0.01*1.05, lambda_M(1)*1.25,text3)
legend([h1 h2], 'Original hazard curve', 'Fitted power-law hazard curve', 'location', 'southwest')
xlabel(Axis_label)
ylabel('Annual rate of exceedance, \lambda')
axis([0.01 figureAxisLimits(2:4)]) % include lower IM values
set(gca, 'xtick', figureXTickVals)




