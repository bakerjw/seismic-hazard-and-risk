% make illustrative loss ratio predictions using HAZUS

% compute loss predictions using HAZUS MH 2.1
% Jack Baker
% 12/20/2017
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
    


pgaVals = 0.01:0.01:1.5; % PGA values of interest

% specify a text case (see fn_HAZUS_loss for allowable options)
analysisCase.codeLevel = 1;
analysisCase.buildingType = 'C2L';
analysisCase.occType = 'COM1';

% perform calculations
[lossRatio, caseLabel, structLossRatio, nonStructAccLossRatio, nonStructDriftLossRatio] = fn_HAZUS_loss(analysisCase, pgaVals);

%% plot contributions of each type of damage
figure
plot(pgaVals, lossRatio, linespec{1}, 'linewidth', 2, 'color', colorspec{1})
hold on
plot(pgaVals, nonStructAccLossRatio, linespec{2}, 'linewidth', 2, 'color', colorspec{2})
plot(pgaVals, structLossRatio, linespec{3}, 'linewidth', 2, 'color', colorspec{3})
plot(pgaVals, nonStructDriftLossRatio, linespec{4}, 'linewidth', 2, 'color', colorspec{4})
xlabel('Peak Ground Acceleration, PGA [g]')
ylabel('Loss ratio')
legend('Total', 'Nonstructural (Acceleration Sensitive)', 'Structural loss', 'Nonstructural (Drift Sensitive)', 'location', 'Northwest')


%% perform calculations for four code levels
clear lossRatio caseLabel

figure
for i=1:4
    analysisCase.codeLevel = 5-i; % flip the ordering
    [lossRatio(:,i), caseLabel{i}] = fn_HAZUS_loss(analysisCase, pgaVals);
end

% manually make shorter labels
shortLabel{1} = 'Pre-code';
shortLabel{2} = 'Low-code';
shortLabel{3} = 'Moderate-code';
shortLabel{4} = 'High-code';


% plot total loss ratios
plot(pgaVals, lossRatio(:,1), linespec{1}, 'linewidth', 2, 'color', colorspec{1})
hold on
plot(pgaVals, lossRatio(:,2), linespec{2}, 'linewidth', 2, 'color', colorspec{2})
plot(pgaVals, lossRatio(:,3), linespec{3}, 'linewidth', 2, 'color', colorspec{3})
plot(pgaVals, lossRatio(:,4), linespec{4}, 'linewidth', 2, 'color', colorspec{4})

xlabel('Peak Ground Acceleration, PGA [g]')
ylabel('Loss ratio')
legend(shortLabel, 'location', 'southeast')

