% example PSHA calculations for book
%
% Created by Jack Baker in ~2008
% Updated March 30, 2016 for use with book
% Revised 4/19/2017 to use PJS Ch 3 rate numbers
% Revised 5/10/2017 to update disaggregation calcs
% Revised 4/4/2019 to add deterministic hazard calculation
% color figure option added 3/10/2021

clear; close all; clc

if 1==0 % greyscale color scheme
    colorspec{1} = [0 0 0];
    colorspec{2} = [0.4 0.4 0.4];
    colorspec{3} = [0.7 0.7 0.7];
    colorspec{4} = [0 0 0];
    
    
    linespec{1} = '-';
    linespec{2} = '-';
    linespec{3} = '-';
    linespec{4} = '--';

    filename_append = ''; % don't add anything to filename

else % color figures
    colorspec{1} = [56 95 150]/255;
    colorspec{2} = [207 89 33]/255;
    colorspec{3} = [158 184 219]/255;
    colorspec{4} = [231 184 0]/255;
    colorspec{5} = [128 0 0]/255;
    
    linespec{1} = '-';
    linespec{2} = '-';
    linespec{3} = '-';
    linespec{4} = '-';
    linespec{5} = '-';
    
    filename_append = '_color';
end

x = logspace(log10(0.001), log10(2), 100); % IM values to consider 
T = 1; % period of interest
IM_label = 'SA(1 s)';
Axis_label = 'Spectral Acceleration, SA(1 s) [g]';
gmpeFlag = 1; % use BJF97


% seismicity parameters
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

% From Table 3.5, \label{tab:grExample_mMax}, fixed rate of M>5, M_max = 8
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
FormatFigureBook
print('-dpdf', ['../figures/example_hazard_GR_source' filename_append '.pdf']); 


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
FormatSubplotFigureBook

set(gcf, 'PaperSize', [6 3.25]);
set(gcf, 'PaperPosition', [0 0 6 3.25]);
print('-dpdf', ['../figures/example_disagg_GR_source' filename_append '.pdf']); 

mBar = [disagg.Mbar disagg2.Mbar] 

% tabulate output
disagg_table = [M_vals' disagg.example' disagg2.example'];

%% disagg conditional on equalling
figure
subplot(1,2,1)
bar(M_vals, disagg.equal, 1,'FaceColor', colorspec{3})
hold on
plot( disagg.equalMbar*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' = ' num2str(x_example) ' g)']);
axis([5 8 0 0.2])
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(1,2,2)
bar(M_vals, disagg2.equal, 1,'FaceColor', colorspec{3})
hold on
plot( disagg2.equalMbar*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' = ' num2str(x_example2) ' g)']);
axis([5 8 0 0.2])
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
FormatSubplotFigureBook

set(gcf, 'PaperSize', [6 3.25]);
set(gcf, 'PaperPosition', [0 0 6 3.25]);
print('-dpdf', ['../figures/example_disagg_equals_GR_source' filename_append '.pdf']); 

mBarEqual = [disagg.equalMbar disagg2.equalMbar]


%% coarsen the disaggregation bins

% extend the arrays to get an even number of values for disaggregation
M_valsEx = [M_vals 8.1];
disagg.exampleEx = [disagg.example 0];
disagg2.exampleEx = [disagg2.example 0];


for j = 1:length(M_valsEx)/2
    idx = 2*j-1; % index value for aggregation
    M_valsCoarse(j) = mean(M_valsEx(idx:idx+1));
    M_disaggCoarse(j) = sum(disagg.exampleEx(idx:idx+1));
    M_disaggCoarse2(j) = sum(disagg2.exampleEx(idx:idx+1));
end

figure
subplot(1,2,1)
bar(M_valsCoarse, M_disaggCoarse, 1,'FaceColor', colorspec{3})
hold on
plot( sum(M_valsCoarse.*M_disaggCoarse)*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' > ' num2str(x_example) ' g)']);
axis([5 8 0 0.4])
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(1,2,2)
bar(M_valsCoarse, M_disaggCoarse2, 1,'FaceColor', colorspec{3})
hold on
plot( sum(M_valsCoarse.*M_disaggCoarse2)*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' > ' num2str(x_example2) ' g)']);
axis([5 8 0 0.4])
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
FormatSubplotFigureBook
set(gcf, 'PaperSize', [6 3.25]);
set(gcf, 'PaperPosition', [0 0 6 3.25]);
print('-dpdf', ['../figures/example_disagg_coarse_GR_source' filename_append '.pdf']); 

mBarCoarse = [sum(M_valsCoarse.*M_disaggCoarse) sum(M_valsCoarse.*M_disaggCoarse2)]


%% disaggregation with epilson 

figure
h1 = subplot(1,2,1)
bar(M_vals, disagg.M_Eps, 1, 'stacked')
colormap(h1, gray)
hold on
% plot( sum(M_vals.*disagg.equal)*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' > ' num2str(x_example) ' g)']);
axis([5 8 0 0.2])
text(-0.1,-0.07,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(1,2,2)
bar(M_vals, disagg2.M_Eps, 1, 'stacked')
hold on
% plot( sum(M_vals.*disagg2.equal)*[1 1], [0 1], ':k', 'linewidth', 2) % plot mean magnitude
hx = xlabel('Magnitude, M');
hy = ylabel(['P(m | ' IM_label ' > ' num2str(x_example2) ' g)']);
axis([5 8 0 0.2])
c = colorbar;
c.TickLabels = num2str(disagg.epsVals');
c.Label.String = '\epsilon';
colormap gray
text(-0.1,-0.07,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
FormatSubplotFigureBook

set(gcf, 'PaperSize', [6 3.25]);
set(gcf, 'PaperPosition', [0 0 6 3.25]);
print('-dpdf', ['../figures/example_disagg_Epsilon_GR_source' filename_append '.pdf']); 

epsBar = [disagg.epsBar disagg2.epsBar]

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

FormatSubplotFigureBook

set(gcf, 'PaperSize', [6 3.25]);
set(gcf, 'PaperPosition', [0 0 6 3.25]);
print('-dpdf', ['../figures/example_hazard_curve_derivative' filename_append '.pdf']); 


%% summary plot
figure
h1 = loglog(x, lambda.x,'-', 'linewidth', 2, 'color', colorspec{1});
hold on
% plot(x_example, lambda.example, 'ok')
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
% axis(figureAxisLimits)
axis([0.01 figureAxisLimits(2:4)]) % include lower IM values
set(gca, 'xtick', figureXTickVals)
FormatFigureBook
print('-dpdf', ['../figures/example_hazard_metrics' filename_append '.pdf']); 


%% hazard curve plus deterministic amplitude

M_deterministic = max(M_vals);
[sa, sigma] = gmpe_eval(T, M_deterministic, rup, gmpeFlag);
sa84 = exp(log(sa)+sigma)
deterministicRate = exp(interp1(x,log(lambda.x),sa84))


figure
h1 = loglog(x, lambda.x,'-', 'linewidth', 2, 'color', colorspec{1});
hold on
h2 = plot([sa84 sa84], figureAxisLimits(3:4),'-', 'linewidth', 2,'Color', colorspec{3});
h3 = plot([figureAxisLimits(1) sa84], deterministicRate*[1 1], ':k', 'linewidth', 1);
legend([h1 h2], 'Probabilistic hazard curve', ['Deterministic amplitude'], 'location', 'northeast')
xlabel(Axis_label)
ylabel('Annual rate of exceedance, \lambda')
axis(figureAxisLimits) 
set(gca, 'xtick', figureXTickVals)
FormatFigureBook
print('-dpdf', ['../figures/example_DSHA' filename_append '.pdf']); 


