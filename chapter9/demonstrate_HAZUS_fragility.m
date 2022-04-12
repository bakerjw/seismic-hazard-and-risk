% Use HAZUS fragility parameters to plot example fragility functions
% Jack Baker
% modified 1/7/2021 to allow optional color figures

% 
% INPUT VARIABLES:
%
%   analysisCase.codeLevel      % flag for code level:
%                                   1 --> High Code
%                                   2 --> Moderate Code
%                                   3 --> Low Code
%                                   4 --> Pre-Code
%
%   analysisCase.buildingType 	% 2- or 3-letter code for construction type. Allowable options:
%                                     W1 	Wood, Light Frame (< 5,000 sq. ft.) 
%                                     W2 	Wood, Commercial and Industrial (> 5,000 sq. ft.)
%                                     S1L 	Steel Moment Frame 
%                                     S1M 	Steel Moment Frame 
%                                     S1H 	Steel Moment Frame 
%                                     S2L 	Steel Braced Frame 
%                                     S2M 	Steel Braced Frame 
%                                     S2H 	Steel Braced Frame 
%                                     S3 	Steel Light Frame 
%                                     S4L 	Steel Frame with Cast?in?Place Concrete Shear Walls 
%                                     S4M 	Steel Frame with Cast?in?Place Concrete Shear Walls 
%                                     S4H 	Steel Frame with Cast?in?Place Concrete Shear Walls 
%                                     S5L 	Steel Frame with Unreinforced Masonry Infill Walls 
%                                     S5M 	Steel Frame with Unreinforced Masonry Infill Walls 
%                                     S5H 	Steel Frame with Unreinforced Masonry Infill Walls 
%                                     C1L 	Concrete Moment Frame 
%                                     C1M 	Concrete Moment Frame 
%                                     C1H 	Concrete Moment Frame 
%                                     C2L 	Concrete Shear Walls 
%                                     C2M 	Concrete Shear Walls 
%                                     C2H 	Concrete Shear Walls 
%                                     C3L 	Concrete Frame with Unreinforced Masonry Infill Walls 
%                                     C3M 	Concrete Frame with Unreinforced Masonry Infill Walls 
%                                     C3H 	Concrete Frame with Unreinforced Masonry Infill Walls 
%                                     PC1 	Precast Concrete Tilt?Up Walls 
%                                     PC2L 	Precast Concrete Frames with Concrete Shear Walls 
%                                     PC2M 	Precast Concrete Frames with Concrete Shear Walls 
%                                     PC2H 	Precast Concrete Frames with Concrete Shear Walls 
%                                     RM1L 	Reinforced Masonry Bearing Walls with Wood or Metal Deck Diaphragms 
%                                     RM1M 	Reinforced Masonry Bearing Walls with Wood or Metal Deck Diaphragms 
%                                     RM2L 	Reinforced Masonry Bearing Walls with Precast Concrete Diaphragms 
%                                     RM2M 	Reinforced Masonry Bearing Walls with Precast Concrete Diaphragms 
%                                     RM2H 	Reinforced Masonry Bearing Walls with Precast Concrete Diaphragms 
%                                     URML 	Unreinforced Masonry Bearing Walls 
%                                     URMM 	Unreinforced Masonry Bearing Walls 
%                                     MH 	Mobile Homes 
%
%   analysisCase.occType  	% 4- or 5-letter code for construction type. Allowable options:
%                                     RES1	Single Family Dwelling
%                                     RES2	Mobile Home
%                                     RES3	Multi Family Dwelling
%                                     RES4	Temporary Lodging
%                                     RES5	Institutional Dormitory
%                                     RES6	Nursing Home
%                                     COM1	Retail Trade
%                                     COM2	Wholesale Trade
%                                     COM3	Personal and Repair Services
%                                     COM4	Professional/Technical/ Business  Services
%                                     COM5	Banks/Financial Institutions
%                                     COM6	Hospital
%                                     COM7	Medical Office/Clinic
%                                     COM8	Entertainment & Recreation
%                                     COM9	Theaters
%                                     COM10	Parking
%                                     IND1	Heavy
%                                     IND2	Light
%                                     IND3	Food/Drugs/Chemicals
%                                     IND4	Metals/Minerals Processing
%                                     IND5	High Technology
%                                     IND6	Construction
%                                     AGR1	Agriculture
%                                     REL1	Church/Membership/Organization
%                                     GOV1	General Services
%                                     GOV2	Emergency Response
%                                     EDU1	Schools/Libraries
%                                     EDU2	Colleges/Universities
%
%   pgaVals                    PGA values for which to compute loss ratios
%

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
    

% specify a text case (see fn_HAZUS_loss for allowable options)
analysisCase.codeLevel = 1;
analysisCase.buildingType = 'C2L';
analysisCase.occType = 'COM1';
pgaVals = 0.01:0.01:1.5; % PGA values of interest
pgaVals = pgaVals(:); % make sure pgaVals is a column 
pgaEx = 0.5; % example value


load hazusData % load data created by import_HAZUS_data.m

% find index for building type
idxBldg = find(strcmp(analysisCase.buildingType, hazusData.buildingTypeCode));
% find index for occupancy type
idxOcc = find(strcmp(analysisCase.occType, hazusData.occCode));


% get fragility parameters for the given structure type and code level
medianDS = hazusData.medians{analysisCase.codeLevel}(idxBldg,:);
betaDS = hazusData.betas{analysisCase.codeLevel}(idxBldg,:);
assert(~isnan(medianDS(1)), ['Error, this building type and code level is not allowed']) % make sure an appropriate occupancy type was specified

% example numbers
p1 = normcdf(log(pgaEx), log(medianDS(1)), betaDS(1))
p2 = normcdf(log(pgaEx), log(medianDS(2)), betaDS(2))
p1_Equal = p1-p2


figure
hold on
for i = 1:4
    plot(pgaVals, normcdf(log(pgaVals), log(medianDS(i)), betaDS(i)), linespec{i}, 'linewidth', 2, 'color', colorspec{i})
end
plot([pgaEx pgaEx], [0 1], ':k')
plot([0 pgaEx], [p1 p1], ':k')
plot([0 pgaEx], [p2 p2], ':k')
legend('ds_1', 'ds_2', 'ds_3', 'ds_4', 'location', 'Northeast')
xlabel('Peak Ground Acceleration, PGA [g]')
ylabel('P(DS \geq ds_i | PGA = x)')
axis([0 1.5 0 1])

% annotate text results for example cases
text1 = ['P(DS \geq ds_1) = ' num2str(p1,2)];
text2 = ['P(DS \geq ds_2) = ' num2str(p2,2)];
text(.05, p1+0.03,text1,'FontSize',7)
text(.05, p2+0.03,text2,'FontSize',7)
text(.55, p2+0.15,['P(DS = ds_1) = ' num2str(p1_Equal,1)],'FontSize',7)




%%  loss ratios for the given occupancy type
lossStruct  = hazusData.lossStruct(idxOcc,:);
lossAccNS   = hazusData.lossAccNS(idxOcc,:);
lossDriftNS = hazusData.lossDriftNS(idxOcc,:);

% example numbers for text
lossTotal = (lossStruct + lossAccNS + lossDriftNS)/100 % loss ratios given DS

diff(-[1 normcdf(log(pgaEx), log(medianDS), betaDS) 0]) % probabilities of each DS, given PGA


%  calculate loss ratio for all PGAs
[lossRatio] = fn_HAZUS_loss(analysisCase, pgaVals);
lossEx = fn_HAZUS_loss(analysisCase, pgaEx) % value for specific case (without rounding)

% plot loss ratio
figure
% plot(pgaVals, structLossRatio, '--')
% hold on
% plot(pgaVals, nonStructAccLossRatio, ':')
% plot(pgaVals, nonStructDriftLossRatio, '-.')
plot(pgaVals, lossRatio, '-', 'linewidth', 2, 'color', colorspec{1})
hold on
plot([pgaEx pgaEx 0], [0 lossEx lossEx], ':k', 'linewidth', 1)
plot(pgaEx, lossEx, 'ok', 'linewidth', 1)
xlabel('Peak Ground Acceleration, PGA [g]')
ylabel('Mean loss ratio, E[C | PGA]')
text(0.56, 0.12, sprintf('PGA = %.2f g \nLoss Ratio = %.2f', pgaEx, lossEx), 'FontSize', 9)




