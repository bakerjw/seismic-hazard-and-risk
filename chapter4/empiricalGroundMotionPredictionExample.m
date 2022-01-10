%Evaluate the BJF97 model for several M, R cases and write the results (Figure 4.26, Table 4.1)
%Brendon Bradley
%May 2020

T=1; %Vibration period to consider (s)
if T > 2; disp('Error: max(T) = 2s for the BJF97 model'); return; end

%earthquake scenario(s) to consider
M = 6.5; %moment magnitude
R = [3,10,30]; %source-to-site distances (km)
Sa1pt0_threshold = 0.4; %IM threshold to compute exceedance probabilities for


%get the GMM prediction - note this requires to GMM to be in the same
%directory as this file, or use the 'addpath' function to point to another
%directory. more information on this function can be found at:
%https://github.com/bakerjw/seismic-hazard-and-risk/tree/main/utils/gmm_bjf97.m
Fault_Type=1; 
Vs=500;
for i=1:length(R)
    [sa(i), sigma(i)] = gmpe_bjf97(M, R(i), T, Fault_Type, Vs);
end

%output the numerical values to the command line 
fprintf('Output values \n');
fprintf('--------------------------------------------------------- \n');
fprintf('     R[km]   mu_lnIM   im_50 [g] sigma_lnIM  Prob[SA(1s) > 0.4g] \n');
fprintf('--------------------------------------------------------- \n');
for i=1:length(R)
    fprintf('%10.3f %10.3f %10.3e %10.3f %10.3e \n',R(i),log(sa(i)),sa(i),sigma(i),1-normcdf(log(Sa1pt0_threshold),log(sa(i)),sigma(i)));
end
fprintf('--------------------------------------------------------- \n');

%Plotting preparation
%get for the full R value range to make plot
R_plot=1:1:100;
for j=1:length(R_plot)
    [sa_plot(j), sigma_plot(j)] = gmpe_bjf97(M, R_plot(j), T, Fault_Type, Vs);
end

%get the necessary details for plotting the pdf distribution
eps_range=-3.7:0.01:3.7; %range of epsilon values to plot
for i=1:length(R)
    sa_vals(:,i)=sa(i)*exp(eps_range*sigma(i));
    sa_pdf(:,i)=normpdf(log(sa_vals(:,i)),log(sa(i)),sigma(i));
    sa_pdf_norm(:,i)=0.6*sa_pdf(:,i)/max(sa_pdf(:,i));
end

%details for markers
Colors=[1 0 0;0 0 1;0 1 0]; %
Markers={'o';'sq';'>'};

%make figure
figure;
%plot BJF97 prediction for full R range
h(2)=loglog(R_plot,sa_plot,'-k','LineWidth',2); hold on;
h(3)=loglog(R_plot,sa_plot.*exp(sigma_plot),'--k','LineWidth',2);
loglog(R_plot,sa_plot.*exp(-sigma_plot),'--k','LineWidth',2);
%plot the SA(1s) threshold
h(1)=loglog(R_plot,ones(length(R_plot),1)*Sa1pt0_threshold,':','Color',[0.5 0.5 0.5],'LineWidth',2);
%plot the lognormal distributions
for i=1:length(R)
    loglog(R(i)*ones(length(eps_range),1),sa_vals(:,i),'-','Color',[0.5 0.5 0.5],'LineWidth',1);
    loglog(R(i)*(1+sa_pdf_norm(:,i)),sa_vals(:,i),'-','Color',[0.5 0.5 0.5],'LineWidth',1);
    %fill in the values which exceed the threshold
    for j=1:length(sa_pdf_norm)
        if (sa_vals(j,i)>Sa1pt0_threshold)
            loglog(R(i)*[1 1+sa_pdf_norm(j,i)],sa_vals(j,i)*[1 1],'-','Color',[0.5 0.5 0.5]);
        end
    end
end

%add text
set(groot,'defaultLegendInterpreter','latex')
legend(h,'Target SA(1 s)','$\mu_{\ln SA}$','$\mu_{\ln SA} \pm \sigma_{\ln SA}$','Location','SouthWest');

%axes formatting
xlim([1 100]); 
ylim([0.01 2]);
xlabel('Source-to-site distance, R [km]'); ylabel('Spectral acceleration, SA(1 s) [g]');


