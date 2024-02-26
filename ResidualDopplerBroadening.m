% Residual Doppler Broadening Calculator
% Seth Erickson
% Based on convolving the Gaussian Doppler lineshape with the natural
% Lorentzian lineshape for each pair of comb teeth which builds up the
% excitation

delta=-4E12:1E9:4E12;
nu=-3E6:1E3:3E6;
dD=delta(2)-delta(1);

FWHM=2000E9;
sig=FWHM/(2*sqrt(2*log(2)));
I = @(d) (1/(sig*sqrt(2*pi)))*exp(-d^2/(2*sig^2));

delta=delta(delta~=0);
len=length(delta);

kb=1.38E-23;
M=1.44E-25;
T=373;
c=2.998E8;

a=c*sqrt(M/(2*pi*kb*T));
b=(M*c^2)/(2*kb*T);

S= @(d) (a/abs(d))*exp(-b*nu.^2/d.^2);
gam=330E3;
lor=(1/pi)*(gam./2)./(nu.^2+gam.^2/4);

Sum=0;
j=0;
for i = delta
    Sum=Sum+conv(S(i),lor,'same')*(I(i).^2);
    j=j+1;
    j/len
end
Sum=Sum/len;
Sum=Sum/max(Sum);

figure(1)
plot(nu,Sum,'k')
xlabel('Frequency (Hz)')
ylabel('$S_{total}(\nu)$','interpreter','latex')
set(gca,'fontname','times')
set(gca,'fontsize',24)
set(findall(gcf,'Type','line'),'LineWidth',2)

