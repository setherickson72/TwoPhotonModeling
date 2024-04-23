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
FWHM=75E9; % BW is std, so put in FWHM here
bw=FWHM./(2*sqrt(2*log(2)));
frep=1E8;

freqs=-5*bw:frep:5*bw;
Spectrum=frep*Itot/(bw*sqrt(pi))*exp(-freqs.^2/bw.^2);

rhos=[];
Ls=[];

Es=sqrt(2*Spectrum/(c*eps0));

j=1;
k=1;
jmax=200;
kmax=200;
wb=waitbar((k-1)/kmax,['Completion = ',num2str(100*k/kmax,3),'%']);
gdd=linspace(-4E8*1E-30,4E8*1E-30,kmax);
Z=linspace(-2.5E-3,2.5E-3,jmax);
dz=Z(2)-Z(1);

for GDD=gdd
    chirp=exp(1i*pi^2*GDD.*(freqs.^2));
 
    for z=Z
        
    EsZf=Es.*exp(1i*2*pi*freqs*z/c).*chirp;
    EsZb=Es.*exp(-1i*2*pi*freqs*z/c).*chirp;
    x1Sq=(dgk/hb)^2*sum(abs(EsZf).^2);
    x1x2=(dgk*dke/hb^2)*sum(EsZf.*flip(EsZb));
    x2Sq=(dke/hb)^2*sum(abs(EsZb).^2);
    
    xeff=x1x2/(2*Del);
    del=(x1Sq-x2Sq)/(4*Del);
    
    rhoSS=(1/2)*(abs(xeff)^2)/(abs(xeff)^2+gam^2+4*del^2);
    Ls=[Ls,z];
    rhos(j,k)=rhoSS;
    
    j=j+1;
    end
    k=k+1;
    j=1;
    waitbar((k-1)/kmax,wb,['Completion = ',num2str(100*(k-1)/kmax,3),'%'])
end
close(wb)

figure(3)
imagesc(gdd*1E30,Z*1E3,rhos.*gam)
colormap turbo
colorbar
clim([0,round(max(max(rhos.*gam)),2,'significant')])
xlabel('GDD (fs$^2$)','interpreter','latex')
ylabel('Distance from pulse center (mm)')
set(gca,'FontSize',28,'FontName','Times')

figure(4)
plot(gdd*1E30,trapz(rhos.*gam.*dz,1),'k')
xlabel('GDD (fs$^2$)','interpreter','latex')
ylabel('Total Fluorescence (arb)')
xlim('tight')
set(gca,'FontSize',28,'FontName','Times')
