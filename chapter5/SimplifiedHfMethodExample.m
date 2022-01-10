function SimplifiedHfMethodExample
%Plot several figures (Figure 5.15, .16, .17, .19) which illustrate the
%essential components of the simplified-physics 'stochastic' method.
%Brendon Bradley
%June 2020

%figure to produce
runtype=4; %1=source - Figure 5.15; 
%           2=path - Figure 5.16; 
%           3=site - Figure 5.17;
%           4=convolution to create complete FAS and time series - Figure 5.19


if runtype==1 %source function
    
    %magnitude and stress parameter values
    M=[6.0 6.0 7.5 7.5];
    DSigma=[50 100 50 100]; %bar    
    
    for i=1:length(M)
        dSigma=DSigma(i);
        f=logspace(-2,2,100);        
        E=sourceSpectra(f,M(i),dSigma);
        
    end
    
    %Plotting details
    figure
    %colors and line styles
    Colors=[0.5 0.5 0 0];
    LineStyles={'-';'--';'-';'--'};
    %plot source spectra
    for i=1:length(M)
        loglog(f,E,'LineStyle',LineStyles{i},'Color',Colors(i)*[1 1 1],'LineWidth',2); hold on;
    end
    xlabel('Frequency, f [Hz]'); ylabel('Fourier Acceleration Spectrum [m/s^2/Hz]');
    xlims=[0.01 10]; ylims=[1e-2 2.5e1];
    xlim(xlims); ylim(ylims);
    text(1.3e-2,3e-1,'M7.5','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0]);
    text(1e-1,5e-2,'M6.0','HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0]);
    text(2,3.0,'\Delta\sigma=10 MPa','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',14,'Color',[0 0 0]);
    text(2,1.6,'\Delta\sigma=5 MPa','HorizontalAlignment','left','VerticalAlignment','top','FontSize',14,'Color',[0 0 0]);
    
elseif runtype==2 %path attenuation
    %parameters for calculation
    f=[0.1 5]; %frequencies of interest
    V=3.6; %shear-wave velocity at source location (km/s)
    Q=[180 371]; %anelastic attenuation values considered for the respective frequencies
    region=2; %1=CEUS; 2=WUS
    
    %source-to-site distances to consider
    R=logspace(log10(1.25),log10(1000/1.25),100);
    
    %get geo spreading function
    Z=geospreading(R,region);
    
    %get path attenuation
    for i=1:length(f)
        P(:,i)=Z.*exp(-pi*f(i)*R/(Q(i)*V));
    end
    
    %plotting
    Colors=[0 0.3 0.6];
    LineStyles={'-';'--';'-.'};
    figure;         
    loglog(R,Z,'LineStyle',LineStyles{1},'Color',Colors(i)*[1 1 1],'LineWidth',2); hold on;
    for i=1:length(f)
        loglog(R,P(:,i),'LineStyle',LineStyles{i+1},'Color',Colors(i+1)*[1 1 1],'LineWidth',2); hold on;
    end
    xlabel('Source-to-site distance, R [km]'); ylabel('Path attenuation, P');
    xlims=[1 1e3]; ylims=[1.e-3 1.e0];
    xlim(xlims); ylim(ylims);
    %plot text
    text(35,4e-2,{'Geometric';'spreading'},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0]);
    bodyWaveLineX=[10 15];
    plot(bodyWaveLineX,1./bodyWaveLineX*1.2,'-k');
    text(geomean(bodyWaveLineX),geomean(1./bodyWaveLineX)*1.2,'1/R','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',14,'Color',[0 0 0],'Interpreter','latex');
    surfWaveLineX=[180 180*1.5];
    plot(surfWaveLineX,(1/70)*(130./surfWaveLineX).^0.5*1.2,'-k');
    text(geomean(surfWaveLineX),geomean((1/70)*(130./surfWaveLineX).^0.5)*1.2,'$1/{\surd}R$','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',14,'Color',[0 0 0],'Interpreter','latex');
    text(600,4e-3,{'f=0.1 Hz';'Q=180'},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14,'Color',[0 0 0]);
    text(100,4e-3,{'f=5 Hz';'Q=371'},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',14,'Color',[0 0 0]);
    
elseif runtype==3 %site amplification
    
    %crustal velocity profile
    [z,Vs,rho]=booreSiteAmpProfile;
    %QWL-based amplification
    Vsrc=3.5; rhosrc=2.8; %velocity and density at the source location (km/s and t/m^3, resp.)
    [f,Amp]=QWL(z,Vs,rho,Vsrc,rhosrc); %amplification based on quarter-wavelength theory
    
    %add dimuniation due to kappa
    kappa0=[0 0.01 0.02 0.04 0.08];
    for i=1:length(kappa0)
        G(:,i) = Amp.*exp(-pi*kappa0(i)*f);
    end
        
    %plot the amplification
    figure;
    loglog(f,G,'-k','LineWidth',2); hold on;
    xlabel('Frequency, f [Hz]'); ylabel('Site Amplification, S');
    xlims=[5e-2 1e2]; ylims=[0.5 5];
    xlim(xlims); ylim(ylims);
    YTicks=[0.5:0.1:5];
    YTickLabels=cell(length(YTicks),1);
    YTickLabels{6}='1'; YTickLabels{16}='2';YTickLabels{26}='3';
    YTickLabels{36}='4';YTickLabels{46}='5';
    set(gca,'YTick',YTicks,'YTickLabel',YTickLabels);
    text(20,2.5,'$\kappa_0=0.0$','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0],'Interpreter','latex');
    text(17,1.65,'0.01','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',14,'Color',[0 0 0],'Interpreter','latex');
    text(9,1.5,'0.02','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0],'Interpreter','latex');
    text(5,1.3,'0.04','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0],'Interpreter','latex');
    text(2,1.05,'0.08','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0],'Interpreter','latex');
    
elseif runtype==4 %combined source+path+site for scenario rupture
    
    %step 1: Compute the FAS model
    %source
    M=6;
    dSigma=50;
    f=logspace(-2,2,100);
    E=sourceSpectra(f,M,dSigma);
    %path
    R=20;
    region=2; %1=CEUS; 2=WUS
    G=geospreading(R,region);
    V=3.5;
    for i=1:length(f)
        if f(i)<1
            Q=180;
        else
            Q=180*f(i)^0.45;
        end
    end
    D=exp(-pi*f*R./(Q*V));
    P=G.*D;
    
    %site
    %get profile
    [z,Vs,rho]=booreSiteAmpProfile;
    %get QWL-based amplification
    Vsrc=3.5; rhosrc=2.8;
    [f_amp,Amp]=QWL(z,Vs,rho,Vsrc,rhosrc);
    %interpolate onto f array
    [I]=interpAmpFn(f_amp,Amp,f);
    %add dimuniation due to kappa
    kappa0=0.045;
    K=exp(-pi*kappa0*f);
    S = I.*K;
    
    %total spectrum
    A=E.*P.*S;
    
    %step 2: white noise generation, and windowing
    %use Boore and Thompson (2014) to get time series duration
    fa = 10^(2.181-0.496*M);
    Ds=1/(2*fa);
    Dp=PathDurationBT14(R);
    Tgm=Ds+Dp;
    
    %create the time domain white noise
    dt=0.005;
    t=0:dt:min(3*Tgm,26);
    tshift=3;
    rng(1); %set random seed to be repeatable.
    noise=randn(1,length(t));
    
    %Saragoni and Hart windowing function
    eps=0.2; eta=0.05;
    b= -(eps*log(eta))/(1+eps*(log(eps)-1));
    c = b/eps;
    a = (exp(1)/eps)^b;
    ftgm=2;
    teta=ftgm*Tgm;          
    w=a*(t/teta).^b.*exp(-c.*(t/teta));
    
    %now plot windowed noise
    figure; 
    normCoeff=max(abs(noise));    
    plot(t,w,'LineWidth',2,'Color',0.7*[1 1 1]); hold on;
    plot(t,-w,'LineWidth',2,'Color',0.7*[1 1 1]);
    plot(t,noise.*w/normCoeff,'k','LineWidth',1.); hold on
    plot([-1.0 0],[0 0],'-k');
%     xlabel('Time, t [s]'); 
    ylabel('Acceleration (unscaled)');
    text(8,0.5,{'Windowing function'},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',0.7*[1 1 1]);
    text(24.5,0.85,{'Step 2: Windowed white noise'},'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0]);
    xlims=[-1.5 25]; ylims=[-1.05 1.05];
    xlim(xlims); ylim(ylims);
    set(gca,'XTickLabel',{});
    
    %step 3: shape the windowed noise and convert back to time domain
    %first convert into fourier spectrum
    [f_,fas]=getFourierSpectra(noise.*w,dt);
    norm=sqrt(sum(real(fas).^2)/length(fas));
    fasN=fas/norm;
    
    %now the inverse fft to get the time series
    nfft=2^ceil(log(length(noise.*w))/log(2));
    Fas=fft(noise.*w,nfft);
    Fas_orig=Fas;
    rms_norm=sqrt(sum(real(Fas).^2)/length(Fas));
    df=(1/dt)*(1)/nfft;
    for i=2:nfft/2
        f_(i)=(i-1)*df;
        [A_(i)]=interpFourierSpectrum(f,A,f_(i));
        sfact(i)=A_(i)/rms_norm;
        Fas(i)=sfact(i)*Fas(i);
        Fas(nfft+1-i)=sfact(i)*Fas(nfft-i);
    end
    sfact(1)=0;
    Fas(1)=0;
    fNq=1/(2*dt);
    [ANq]=interpFourierSpectrum(f,A,fNq);
    Fas(nfft)=ANq*Fas(nfft)/norm;
    %plot
    figure
    loglog(f_(2:end),abs(Fas(2:nfft/2)),'k','LineWidth',1.); hold on;
    loglog(f,A,'LineWidth',2,'Color',0.7*[1 1 1]); hold on;    
    xlabel('Frequency, f [Hz]'); ylabel('Fourier Acceleration Spectrum [m/s^2/Hz]');
    xlims=[0.01 100]; ylims=[5e-6 5e-1];
    xlim(xlims); ylim(ylims);
    text(xlims(1)*1.2,ylims(2)*0.7,{'Step 4: Spectrum-scaled noise'},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0]);
    text(0.04,1e-3,{'Modelled spectrum, A(f)'},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14,'Color',0.7*[1 1 1]);
    
    %comparison plot;
%     figure;
%     plot([1:nfft/2],sfact);
%     semilogy([1:nfft],abs(Fas_orig),[1:nfft],abs(Fas)); hold on;
    
%     figure; plot(f_,A_);
%     Atotal=[A_ flipud(A_)];
%     figure; plot([f_ (200-f_)],Atotal);
    acc_scaled=(1/dt)*ifft(Fas,nfft);
    t2=[0:1:nfft-1]*dt;
    figure;
    plot(t2(1:length(t)),acc_scaled(1:length(t)),'-k','LineWidth',1.); hold on; 
    plot([-1 0],[0 0],'-k','LineWidth',1.);
    xlabel('Time, t [s]'); ylabel('Acceleration, a [m/s^2]');
    xlims=[-1.5 25]; ylims=[-0.7 0.7];
    xlim(xlims); ylim(ylims);
    text(24.5,0.6,{'Step 5: Simulated ground motion'},'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14,'Color',[0 0 0]);
    
    
end

end

function [A_target]=interpFourierSpectrum(f,A,f_target)
%interpolate with constants at bounds

for i=1:length(f_target)
   if f_target(i)<f(1)
       A_target(i)=A(1);
   elseif f_target(i)>f(end)
       A_target(i)=A(end);
   else
       A_target(i)=interp1(f,A,f_target(i));
   end
end

end

function [f,fas]=getFourierSpectra(acc,dt)

nfft=2^ceil(log(length(acc))/log(2));
Fas=fft(acc,nfft);
% figure;plot([1:1:nfft],abs(Fas));
f=(1/dt)*(0:nfft/2)/nfft;
fas=Fas(1:length(f));
%now remove the first point for f=0
f(1)=[];
fas(1)=[];

%testing;
% acc2=ifft(Fas);
% A=1;

end

function Dp=PathDurationBT14(R)

tabulated=[
0 0
7 2.4
45 8.4
125 10.9
175 17.4
270 34.2
];

for i=1:length(R)
   if R(i)<270
        Dp(i)=interp1(tabulated(:,1),tabulated(:,2),R(i));
   else
       Dp(i)=34.2+0.156*(R(i)-270);
   end
end

end



function [A]=interpAmpFn(f_amp,Amp,f)
%interpolate with constants at bounds

for i=1:length(f)
   if f(i)>f_amp(2)
       A(i)=Amp(2);
   elseif f(i)<f_amp(end)
       A(i)=Amp(end);
   else
       A(i)=interp1(f_amp(2:end),Amp(2:end),f(i));
   end
end

end


function E=sourceSpectra(f,M,dSigma)

betaS=3.5; rhoS=2.8; %units: km/s and t/m^3
M0=10^((M+6.03)*(3/2)); %N.m
dSigma_2=dSigma*10^5; %stress parameter in Pa
betaS_2=betaS*10^3; %beta in m/s
%corner frequency
f0=(4.9*10^-1)*betaS_2*(dSigma_2/M0)^(1/3);

S=1./(1+(f/f0).^2);

Rad=0.55; %radiation factor
V=1/sqrt(2); %horizontal component partitioning
F=2; %free surface factor
R0=1; %reference distance (km)
%all of the parameters have SI units
%rho in kg/m^3; beta in m/s; R0 in m; 
C=Rad*V*F/(4*pi*(rhoS*1e3)*(betaS*1e3)^3*(R0*1e3));

E=C*M0*S.*(2*pi*f).^2;

end

function Z=geospreading(R,region)

if region==1 %CEUS (Boore 2003 PAGEOPH)
    for i=1:length(R)
        if R(i)<70
            Z(i)=1/R(i);
        elseif R(i)<130
            Z(i)=1/70;
        else
            Z(i)=(1/70)*(130/R(i))^0.5;
        end
    end    
elseif region==2 %WUS (Atkinson + Silva 2000 BSSA)
    for i=1:length(R)
        if R(i)<40
            Z(i)=1/R(i);
        else
            Z(i)=(1/40)*(40/R(i))^0.5;
        end
    end
end

end

function [z,Vs,rho]=booreSiteAmpProfile
%based on Boore (2016, BSSA)
data=[%Z (km)	Vs (km/s)	rho (t/m^3)
0	0.314	1.934
0.0009	0.314	1.934
0.001	0.427	1.989
0.002	0.512	2.024
0.003	0.569	2.046
0.005	0.649	2.075
0.008	0.731	2.103
0.011	0.793	2.122
0.014	0.843	2.136
0.018	0.898	2.152
0.022	0.944	2.163
0.03	1.02	2.184
0.044	1.176	2.222
0.064	1.348	2.26
0.082	1.474	2.287
0.102	1.594	2.312
0.126	1.718	2.338
0.15	1.826	2.36
0.19	1.984	2.392
0.25	2.08	2.412
0.3	2.147	2.426
0.45	2.308	2.459
0.55	2.393	2.477
0.65	2.467	2.493
0.8	2.561	2.513
0.9	2.614	2.524
1	2.663	2.535
1.45	2.839	2.573
2.05	3.014	2.612
2.4	3.094	2.629
2.85	3.18	2.648
3.4	3.271	2.668
4	3.357	2.687
5.2	3.418	2.701
6.65	3.478	2.714
7.85	3.519	2.723
15   3.519 2.723
];

z=data(:,1);
Vs=data(:,2);
rho=data(:,3);

end

function [f,Amp]=QWL(z,Vs,rho,Vsrc,rhosrc)
%applies quarter-wavelength approach to approximate the amplification as a
%function of frequency

%loop over the depths and for each compute the average velocity 
for i=1:length(z)
    if i==1
        Vsbar(1) = Vs(1); rhobar=rho(1);
    else
        dz=diff(z(1:i));
        Vs_avg=Vs(2:i)+diff(Vs(1:i));
        rho_avg=rho(2:i)+diff(rho(1:i));
        Vsbar(i) = z(i) / sum(dz./Vs_avg);
        rhobar(i) = (1/z(i))*sum(rho_avg.*dz);
    end  
end
Zbar = rhobar.*Vsbar;
f=(1/4)*Vsbar./z';

Zsrc = rhosrc*Vsrc;

Amp = sqrt(Zsrc./Zbar);
end
