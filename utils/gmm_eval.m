function [sa, sigma] = gmm_eval(T, M, rup, gmpeFlag)
% master function to call a relevant GMPE and get median Sa and log
% standard deviation
%
% INPUTS
%
% T         IM period of interest
% M         rupture magnitude
% rup           data structure with rupture parameters
% gmpeFlag      =1 for BJF97, =2 for CY14
%
% sa        median spectral acceleration, given rup
% sigma     log standard deviation, given rup
%

if gmpeFlag == 1 % BJF 1997 model
   [sa, sigma] = gmm_bjf97(M, rup.R, T, rup.Fault_Type, rup.Vs30);
    
elseif gmpeFlag == 2 % CY 2014 model
    [sa, sigma] = gmm_cy2014( M, T, rup.R, rup.R, rup.R, rup.Ztor, rup.delta, rup.rupLambda, rup.Z10, rup.Vs30, rup.Fhw, rup.FVS30, rup.region);

else
    fprintf('Invalid gmpeFlag \n')
end


end

