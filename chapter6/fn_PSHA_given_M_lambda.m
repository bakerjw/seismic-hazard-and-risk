function [lambda, example_output, disagg] = fn_PSHA_given_M_lambda(lambda_M, M_vals, T, x, x_example, rup, gmpeFlag)

% Compute PSHA, with rupture rates for each M precomputed

% Created by Jack Baker

%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lambda_M      exceedance rate of EQs for each M
% M_vals        values of M corresponding to lambda_M
% x             IM values of interest
% x_example     example IM value for table
% rup           data structure with rupture parameters
% gmpeFlag      =1 for BJF97, =2 for CY14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




lambda_occur = [-diff(lambda_M) lambda_M(end)]; % find occurence rates from exceedance rates

% p(exceeding each x threshold value | M)
for j = 1:length(x)
    for i = 1:length(M_vals)
        [sa, sigma] = gmpe_eval(T, M_vals(i), rup, gmpeFlag);
        p_given_M(i) = 1 - normcdf(log(x(j)),log(sa),sigma);
    end
    
    lambda.x(j) = sum(lambda_occur .* p_given_M); % rate of exceeding x
    disagg.all(j,:) = (lambda_occur .* p_given_M) / lambda.x(j);
    
end


% calcs for example IM case
for i = 1:length(M_vals)    
    [sa, sigma] = gmpe_eval(T, M_vals(i), rup, gmpeFlag);
    p_ex(i) = 1 - normcdf(log(x_example),log(sa),sigma);
end

example_output = [[1:length(M_vals)]' M_vals' lambda_occur' p_ex' (lambda_occur' .* p_ex')];
lambda.example = sum(lambda_occur .* p_ex);
disagg.example = (lambda_occur .* p_ex) / lambda.example;
disagg.Mbar = sum(M_vals.*disagg.example);

% disagg conditional on occurence for example IM case
xInc = x_example*1.02; % do computations at an increment on x
for i = 1:length(M_vals)    
    [sa, sigma] = gmpe_eval(T, M_vals(i), rup, gmpeFlag);
    pInc(i) = 1 - normcdf(log(xInc),log(sa),sigma);
end
lambdaInc = sum(lambda_occur .* pInc);
disagg.equal = ((lambda_occur .* p_ex) - (lambda_occur .* pInc)) / (lambda.example-lambdaInc);
disagg.equalMbar = sum(M_vals.*disagg.equal);


% disaggs with epsilon
deltaEps = 1; % final binning
epsVals = -3:deltaEps:3; % epsilon bins

deltaEpsFine = 0.1; % initial finer binning
epsValsFine = -3.5:deltaEpsFine:3.5; % midpoints of bins
p_eps = normpdf(epsValsFine) * deltaEpsFine; % estimate PDF using a PMF with discrete epsilon increments
lambda_M_and_eps = lambda_occur' * p_eps; % rate of events with a given magnitude and epsilon  


for i = 1:length(M_vals)
    [sa, sigma] = gmpe_eval(T, M_vals(i), rup, gmpeFlag);    
    Ind(i,:) = (log(sa) + epsValsFine*sigma > log(x_example)); % indicator that the M/epsilon value causes IM > x_example  
end
exceedRatesFine = Ind .* lambda_M_and_eps; % rates of given M/epsilon values exceeding IM
lambdaExceed = sum(sum(exceedRatesFine)); % this is close to lambda.example, but may differ by a few percent due to epsilon discretization

% compute mean epsilon
epsDeagg = sum(exceedRatesFine) ./ sum(sum(exceedRatesFine));
disagg.epsBar = sum(epsValsFine.*epsDeagg);

% aggregate results to coarser epsilon bins
for j=1:length(epsVals) 
    idx = epsValsFine >= (epsVals(j)-deltaEps/2) & epsValsFine < (epsVals(j)+deltaEps/2);
    exceedRates(:,j) = sum(exceedRatesFine(:,idx),2);
end

disagg.epsVals = epsVals; % return bin midpoints
disagg.M_Eps = exceedRates / lambdaExceed; % magnitude and epsilon disaggregation
disagg.eps = sum(exceedRates) / lambdaExceed; %  epsilon disaggregation


disagg.equalMbar = sum(M_vals.*disagg.equal);

