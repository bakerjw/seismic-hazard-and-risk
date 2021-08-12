function [ lossRatio, caseLabel, structLossRatio, nonStructAccLossRatio, nonStructDriftLossRatio ] = fn_HAZUS_loss( analysisCase, pgaVals )
%
% Function to predict losses using HAZUS MH 2.1
%
% Created by Jack Baker
% June 11, 2016
%
% Note--these are the "Equivalent PGA Structural Fragility Curves from
% section 5.4.4 of Hazus, not the more precise fragilities based on
% capacity curves. They are similar to the capacity curve results, but not
% equivalent.
%
% Note--these are the "Equivalent PGA Structural Fragility Curves from
% section 5.4.4 of Hazus, not the more precise fragilities based on
% capacity curves. They are similar to the capacity curve results, but not
% equivalent.
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
% OUTPUT VARIABLES:
%
%   lossRatio                   loss ratio (total loss) for each PGA
%
%   caseLabel                   text label describing the analysis case
%
%   structLossRatio             loss ratio (structural) for each PGA
%
%   nonStructAccLossRatio       loss ratio (nonstructural acceleration 
%                               sensitive) for each PGA
%
%   nonStructDriftLossRatio     loss ratio (nonstructural drift 
%                               sensitive) for each PGA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load hazusData % load data created by import_HAZUS_data.m

% find index for building type
idxBldg = find(strcmp(analysisCase.buildingType, hazusData.buildingTypeCode));
assert(length(idxBldg)==1, ['Error with specified buiding type, ' analysisCase.buildingType]) % make sure an appropriate building type was specified

% find index for occupancy type
idxOcc = find(strcmp(analysisCase.occType, hazusData.occCode));
assert(length(idxOcc)==1, ['Error with specified occupancy type, ' analysisCase.occType]) % make sure an appropriate occupancy type was specified

% make sure pgaVals is a column 
pgaVals = pgaVals(:);

% get fragility parameters for the given structure type and code level
medianDS = hazusData.medians{analysisCase.codeLevel}(idxBldg,:);
betaDS = hazusData.betas{analysisCase.codeLevel}(idxBldg,:);
assert(~isnan(medianDS(1)), ['Error, this building type and code level is not allowed']) % make sure an appropriate occupancy type was specified

% get loss ratios for the given occupancy type
lossStruct  = hazusData.lossStruct(idxOcc,:);
lossAccNS   = hazusData.lossAccNS(idxOcc,:);
lossDriftNS = hazusData.lossDriftNS(idxOcc,:);

% damage state exceedance probabilities per pgaVals 
for i = 1:length(pgaVals)
    pDsExceed(i,:) = normcdf( log(pgaVals(i)./medianDS)./betaDS);
end
pDsExceed = [pDsExceed zeros(size(pgaVals))]; % pad with zeros for probability of exceeding DS 5, to help with differentiation in the next step

for j = 1:4
    pDsEqual(:,j) = pDsExceed(:,j) - pDsExceed(:,j+1);
end

% compute loss ratios per PGA
structLossRatio         = sum(pDsEqual .* (ones(size(pgaVals)) * lossStruct ),2)/100;
nonStructAccLossRatio   = sum(pDsEqual .* (ones(size(pgaVals)) * lossAccNS  ),2)/100;
nonStructDriftLossRatio = sum(pDsEqual .* (ones(size(pgaVals)) * lossDriftNS),2)/100;

lossRatio = structLossRatio + nonStructAccLossRatio + nonStructDriftLossRatio;

% make a label for the analysis case
caseLabel = [analysisCase.buildingType ', ' analysisCase.occType ', ' hazusData.codeLevel{analysisCase.codeLevel} '-code'];

end

