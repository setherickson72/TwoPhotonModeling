% Steady State Excitation with Comb
% by Seth Erickson

% Constants
c=299792458;
hb=1.05457E-34;
eps0=8.854E-12;

% Constants for our transition
dgk=5.978*8.478E-30; % -e<k|r|g> for our transition
dke=1.98*8.478E-30; % -e<e|r|k> for our transition
Del=2*pi*1E12; % Detuning of the virtual state from real intermediate state
gam=4.193E6; % Excited state decay rate (s-1)

Itot=1E3; %Intensity in Watts/m^2 (1E3=1mW/mm^2)
minbw=10E9; % minimum bandwidth to simulate
maxbw=200E9; % maximum bandwidth to simulate
frep=1E8; % Comb repetition rate

rhos=[];
Ls=[];

freqs=-5*maxbw:frep:5*maxbw;
j=1;
jmax=100; % number of steps to integrate different bandwidths
k=1;
kmax=100; % number of different positions to simulate at each bandwidth

% Added Dispersion
%D2=1E7*1E-30;%23000E-30*5E3;
D2=0;
chirp=exp(1i*pi^2*D2.*(freqs.^2));

wb=waitbar(0,'Progress');
for bw=linspace(minbw,maxbw,kmax)
    Spectrum=frep*Itot/(bw*sqrt(2*pi))*exp(-freqs.^2/(2*bw.^2));
    Es=sqrt(2*Spectrum/(c*eps0));
    Es=Es.*chirp;
    
    for L=linspace(-2.5E-3,2.5E-3,jmax)
    
        EsZf=Es.*exp(1i*2*pi*freqs*L/c);
        EsZb=Es.*exp(-1i*2*pi*freqs*L/c);
        x1Sq=(dgk/hb)^2*sum(abs(EsZf).^2);
        x1x2=(dgk*dke/hb^2)*sum(EsZf.*flip(EsZb));
        x2Sq=(dke/hb)^2*sum(abs(EsZb).^2);
        
        xeff=x1x2/(2*Del);
        del=(x1Sq-x2Sq)/(4*Del);
        
        rhoSS=(1/2)*(abs(xeff)^2)/(abs(xeff)^2+gam^2+4*del^2);
        %Ls=[Ls,L];
        rhos(j,k)=rhoSS;
        j=j+1;
    end
    waitbar(k/kmax,wb,['Progress = ',num2str(100*k/kmax,3),'%'])
    j=1;
    k=k+1;
end
close(wb)

% Plotting
figure(121)
imagesc(linspace(minbw,maxbw,kmax)*1E-9,linspace(-2.5E-3,2.5E-3,jmax)*1E3,rhos*gam)
cb=colorbar;
colormap('jet')
xlabel('Pulse Bandwidth (GHz)')
ylabel('Displacement z (mm)')
set(gca,'FontName','Times')
set(gca,'FontSize',22)
ylabel(cb,'Fluorescence','FontSize',22,'FontName','Times')