% Logic tree calculation for failure rate
% Jack Baker
% 8/3/2018
% modified 1/7/2021 to allow optional color figures

clear; close all; clc;

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
    

%% load hazard curve data from Chapter 6 (ch6/calculations/PSHA_calc_w_epistemic.m)
load hazardEpistemicData
wtHaz = wt; % rename weight vector to avoid ambiguity
clear wt
x = logspace(log10(0.05), log10(2), 100); % IM values to consider 



%% fragility function specification
betaIM = 0.4;
thetaIM = [0.4 0.6];
wtTheta = [0.5 0.5];

for j = 1:length(thetaIM)
    fragility(j,:) = normcdf(log(x), log(thetaIM(j)), betaIM);
end


%% risk calculations over all logic tree branches

idx = 0; % running index for the outputs
for i=1:length(wtHaz) % loop over hazard curve branches
    dLambda =  abs([diff(lambda_x(i,:)) 0]);  % derivative of hazard curve
    
    for j=1:length(thetaIM) % loop over fragility branches
        idx = idx + 1; % increment index
        
        failRate(idx) = sum(fragility(j,:) .* dLambda);
        wtMaster(idx,:) = wtHaz(i) * wtTheta(j);
    end
end

%% mean hazard and mean fragility
for k = 1:length(x)
    lambdaImMean(:,k) = sum(lambda_x(:,k) .* wtHaz); % mean hazard
    fragilityMean(:,k) = sum(fragility(:,k) .* wtTheta'); % mean fragility
end

dLambdaMean =  abs([diff(lambdaImMean) 0]);  % derivative of mean hazard curve
failRateMeanInputs = sum(fragilityMean .* dLambdaMean); % failure rate using mean hazard and mean fragility



%% fractiles of failure rate
% fractiles of hazard
[failRateSort, dataIDX] = sort(failRate); % order branches from lowest rate to highest
weightCum = cumsum(wtMaster(dataIDX)); % order weights appropriately, and cumulatively sum


figure
stairs(failRateSort, weightCum, 'linewidth', 2, 'color', colorspec{3})
hold on
plot(failRateMeanInputs*[1 1], [0 1], '-', 'linewidth', 2, 'color', [0 0 0])
legend('Empirical CDF from logic tree', '\lambda(F) from mean inputs', 'location', 'southeast')
set(gca, 'ylim', [0 1])
xlabel('Annual failure rate, \lambda(F)')
ylabel('Cumulative probability')

%% bar chart of failure rate

xInt = 0.2e-3; % width of intervals to plot
xPlot = xInt/2:xInt:13*xInt; % IM intervals to plot
for i=1:length(xPlot)
    idx = find( (failRate>=xPlot(i)-xInt/2) & (failRate<xPlot(i)+xInt/2));
    yPlot(i) = sum(wtMaster(idx)); % sum weights of branches that fall in this bin
end

figure
h1 = bar(xPlot, yPlot, 1, 'FaceColor', colorspec{3});
hold on
h2 = plot(failRateMeanInputs*[1 1], [0 0.3], '-', 'linewidth', 2, 'color', [0 0 0]);
legend([h1 h2], 'Histogram from logic tree', '\lambda(F) from mean inputs')
xlabel('Annual failure rate, \lambda(F)')
ylabel('Probability')
axis([0 2.5e-3 0 0.3])

%% double-lognormal fragility model

mu_lnTheta = (log(thetaIM(1)) + log(thetaIM(2)))/2
thetaMedian = exp(mu_lnTheta)
betaTheta = sqrt(( (log(thetaIM(1))-mu_lnTheta)^2 + (log((thetaIM(2)))-mu_lnTheta)^2)/2)
betaTot = sqrt(betaIM^2 + betaTheta^2)

failRateDoubleLog = sum(normcdf(log(x), log(thetaMedian), betaTot) .* dLambdaMean) % failure rate using mean hazard and double-lognormal fragility

failRateMedianFragility = sum(normcdf(log(x), log(thetaMedian), betaIM) .* dLambdaMean) % failure rate using mean hazard and fragility with no epistemic

failRateMedianFragility/failRateDoubleLog

%% fragility functions plot
figure
h1 = plot(x, fragility, 'linewidth', 1, 'color', colorspec{2});
hold on
h2 = plot(x,fragilityMean, 'linewidth', 2, 'color', colorspec{3});
h3 = plot(x,normcdf(log(x), log(thetaMedian), betaTot), '--', 'linewidth', 2, 'color', colorspec{1});
legend([h1(1), h2, h3], 'Logic tree branches', 'Mean of logic tree fragilities', 'Double-lognormal model', 'location', 'southeast')
xlabel('Spectral Acceleration, SA(1 s) [g]')
ylabel('P(F | SA=x)')
axis([0 1 0 1])



