% Illustrate risk-targeted hazard level

% Jack Baker
% 12/20/2017, modified 8/1/2018
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
    

%% load previous hazard curve 
x = logspace(log10(0.001), log10(2), 100); % IM values to consider 

% mean hazard from Chapter 6 simplified logic tree 
lambdaMean = [0.0251694954594308,0.0243005973346253,0.0234365472558623,0.0225792462468511,0.0217305158448367,0.0208920837822666,0.0200655715142511,0.0192524837495085,0.0184542000865126,0.0176719688005284,0.0169069027730509,0.0161599775046000,0.0154320311063565,0.0147237661269249,0.0140357530384293,0.0133684351816874,0.0127221349535275,0.0120970610102702,0.0114933162595579,0.0109109064174252,0.0103497489179256,0.00980968197780037,0.00929047363756543,0.00879183062195665,0.00831340688589137,0.00785481173602961,0.00741561744179956,0.00699536627266170,0.00659357691982798,0.00620975028017252,0.00584337459734840,0.00549392996997187,0.00516089224907525,0.00484373635688936,0.00454193906650127,0.00425498128721521,0.00398234990373885,0.00372353921886845,0.00347805204941328,0.00324540052394014,0.00302510662878499,0.00281670254590154,0.00261973082270939,0.00243374441035183,0.00225830660283430,0.00209299090552312,0.00193738085754465,0.00179106982882404,0.00165366080890037,0.00152476620129790,0.00140400763414932,0.00129101579497441,0.00118543029502282,0.00108689956639214,0.000995080793225799,0.000909639876668399,0.000830251431896055,0.000756598814428752,0.000688374172055584,0.000625278518042679,0.000567021820830418,0.000513323105142555,0.000463910559307276,0.000418521643610899,0.000376903194651442,0.000338811520914285,0.000304012485139082,0.000272281569469359,0.000243403919858480,0.000217174366732617,0.000193397419468867,0.000171887232821039,0.000152467544004570,0.000134971579723542,0.000119241932976582,0.000105130410005028,9.24978482376873e-05,8.12139065351012e-05,7.11568294363794e-05,6.22131874592354e-05,5.42775957953396e-05,4.72524139769621e-05,4.10474292661717e-05,3.55795266353024e-05,3.07723482686310e-05,2.65559455227739e-05,2.28664262410104e-05,1.96456002286791e-05,1.68406255681954e-05,1.44036582884222e-05,1.22915077099856e-05,1.04652995717140e-05,8.89014880956513e-06,7.53484361431181e-06,6.37154214337735e-06,5.37548301094735e-06,4.52471043365022e-06,3.79981467097529e-06,3.18368817390031e-06,2.66129764494389e-06];
lambdaIM = lambdaMean; % pick a hazard curve to use

dLambda =  abs([diff(lambdaIM) 0]);  % derivative of hazard curve

%% find risk-targeted fragility

lambdaTarg = 2e-4; % target collapse risk
pTarg = 1-exp(-lambdaTarg*50); % collapse probability in 50 years

betaIM = 0.6; % assumed fragility beta
thetaIM = 0.1:0.01:0.6;
for i=1:length(thetaIM)
    pFail = normcdf(log(x), log(thetaIM(i)), betaIM);
    failRate(i) = sum(pFail .* dLambda);
end
thetaTarget = interp1(failRate,thetaIM,lambdaTarg); % find risk-targeted theta

figure
plot(thetaIM, failRate, '-', 'linewidth', 2, 'color', colorspec{1})
hold on
plot(thetaIM, lambdaTarg*ones(size(thetaIM)), ':k', 'linewidth', 1)
plot(thetaTarget*[1 1], [0 lambdaTarg], '--k', 'linewidth', 1)
legend('\lambda(F | \theta)', 'Target \lambda(F)')
xlabel('Fragility median, \theta [g]')
ylabel('Annual failure rate, \lambda(F)')

% alternate beta
betaIM2 = 0.7; % assumed fragility beta
for i=1:length(thetaIM)
    pFail = normcdf(log(x), log(thetaIM(i)), betaIM2);
    failRate2(i) = sum(pFail .* dLambda);
end
thetaTarget2 = interp1(failRate2,thetaIM,lambdaTarg); % find risk-targeted theta

% another alternate beta
betaIM3 = 0.8; % assumed fragility beta
for i=1:length(thetaIM)
    pFail = normcdf(log(x), log(thetaIM(i)), betaIM3);
    failRate3(i) = sum(pFail .* dLambda);
end
thetaTarget3 = interp1(failRate3,thetaIM,lambdaTarg); % find risk-targeted theta


fig = figure;
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
yyaxis left
plot(x, normcdf(log(x), log(thetaTarget), betaIM), '-', 'linewidth', 2, 'color', colorspec{1})
hold on
plot(x, normcdf(log(x), log(thetaTarget2), betaIM2), '-', 'linewidth', 2, 'color', colorspec{2})
plot(x, normcdf(log(x), log(thetaTarget3), betaIM3), '-', 'linewidth', 2, 'color', colorspec{3})
ylabel('Probability of failure, P(F | SA=x)')
xlabel('Spectral Acceleration, SA(1 s) [g]')
yyaxis right
plot(x, lambdaIM, '--', 'linewidth', 2, 'color', colorspec{4})
ylim([0 0.002])
ylabel('Annual exceedance rate, \lambda(SA>x)')
legend(['\beta = ' num2str(betaIM)], ['\beta = ' num2str(betaIM2)], ['\beta = ' num2str(betaIM3)], 'Ground-motion hazard', 'location', 'southeast')
xlim([0 1])



