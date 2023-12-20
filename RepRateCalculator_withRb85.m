%% Rep Rate Calculator
%   What rep rate will stop your features from overlapping?
%   Seth Erickson

%% User Input
range = [150E6,400E6]; % Range of repetition rates for confusing plot
RepRate=200E6; % Particular Repetition Rate for Spectrum plot (inside of range)
gamma=1000E3; % How broad to model the spectral lines for Spectrum plot

%% Calculations and confusing plot
fClock=385.2845663663E12;
fs=[385.2880175669E12;385.2880096002E12;385.2879981122E12;385.2846002252E12;385.2845922552E12;385.284580778E12]; % Optical frequency for two photon transitions
as=[0.1125,0.146,0.117,0.0125,0.0625,0.175]; % Relative amplitude for each of these transitions

fs85=[385.285142367E12;385.285147084E12;385.2851515948E12;385.2851553976E12;385.285158138E12;385.286664961E12;385.286669461E12;385.286673254E12;385.286675998E12;385.286677451E12];
as85=[0.306,0.167,0.078,0.028,0.006,0.083,0.117,0.111,0.078,0.028];
fs=[fs;fs85];
as=[as*0.27,as85*0.72];

range=range./2; % Needed to account for resonances at each half repRate

n=int32(fClock/range(2)):int32(fClock/range(1));
n=flip(double(n)); % Reverse n so that RepRate is in order
repRates=fClock./n; % Repetition rates where the clock is resonant
RepRates=repRates.*2; % The factor of two comes back to give real RepRates

ns=round(fs*(1./repRates)); % Which n will be closest to any given repRate
fRs=repmat(fs,[1,numel(n)])./ns; % What repRate will be resonant with F=1 -> F'=1
dFs=(fRs-repmat(repRates,[numel(fs),1])).*fClock./repmat(repRates,[numel(fs),1]); % Effective optical frequency seperation

figure(1) % Confusing plot
hold on
for i=1:numel(fs)
    plot(10^-6*RepRates,10^-6.*dFs(i,:),'.')
    i
end
plot(10^-6*RepRates,2*ones(numel(repRates),1),'b--',10^-6*RepRates,-2*ones(numel(repRates),1),'b--')
hold off
xlabel('Repetition Rate (MHz)')
ylabel('\Deltaf_{RR} to next peak')
title('Calculated Repetition Rates')
legend('F=1\rightarrow1','F=1\rightarrow2','F=1\rightarrow3','F=2\rightarrow1','F=2\rightarrow2','F=2\rightarrow3','\pm2 MHz Seperation','')
set(findobj(gcf,'Type','Line'),'LineWidth',1)
set(gca,'FontName','Times')
set(gca,'FontSize',22)

%% Calculate a DCS-ORAFS Spectrum at a Given Repetition Rate

lorentz= @(f,f0,a) (a*(gamma/2)./((f-f0).^2+(gamma/2)^2)); % Define a lorentzian function

fAxis=linspace(-100E6,100E6,10000);

[M,minIn]=min(abs(RepRates-RepRate)); % Identify the index for a given RepRate
dcs=lorentz(fAxis,0,0.375*0.27); % Start spectrum with the clock transition

figure(2)
plot([0,0],[0,1],'--')
hold on
for i = 1:numel(fs)
    plot([dFs(i,minIn),dFs(i,minIn)].*10^-6,[0,1],'--') % Colored lines to mark the transitions
    dcs=dcs+lorentz(fAxis,dFs(i,minIn),as(i)); % Adding florescence from each transition
    i
end

plot(fAxis./10^6,dcs,'k')
xlabel('Relative Optical Frequency (MHz)')
ylabel('Florescence')
xlim([-75,75])
ylim([0,1*10^-6])
legend('F=2\rightarrow4 (clock)','F=1\rightarrow1','F=1\rightarrow2','F=1\rightarrow3','F=2\rightarrow1','F=2\rightarrow2','F=2\rightarrow3','DCS-ORAFS')
set(findobj(gcf,'Type','Line'),'LineWidth',1)
set(gca,'FontName','Times')
set(gca,'FontSize',22)
title(['Repetition Rate = ',num2str(RepRate*10^-6),' MHz'])
subtitle(['Modelled with \Gamma = ',num2str(gamma*10^-6), ' MHz'])

hold off