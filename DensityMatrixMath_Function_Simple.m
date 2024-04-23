% Density Matrix Math
% Seth Erickson

c=299792458;
eps0=8.854E-12;
hb=1.054E-34;
I=1E7;
E=sqrt(2*I/(c*eps0));
dgk=5.978*8.478E-30;
dke=1.98*8.478E-30;
xgk=dgk*E/hb;
xke=dke*E/hb;

a=xgk/2;
c=xke/2;
Del=2*pi*1E12;
del=0;
gam=4.193E6;
i=1i;

dxdt=@(t,X) [-i*a*conj(X(2))+i*conj(a)*X(2)+gam*X(6); ...
     i*a*X(1)+i*Del*X(2)-i*a*X(4)+i*conj(c)*X(3); ...
     i*c*X(2)+i*del*X(3)-i*a*X(5)-gam/2*X(3);...
     -i*a*X(2)+i*conj(a)*conj(X(2))+i*conj(c)*X(5)-i*c*conj(X(5));...
     -i*conj(a)*X(3)+i*c*X(4)+i*(del-Del)*X(5)-i*c*X(6);...
     -i*conj(c)*X(5)+i*c*conj(X(5))-gam*X(6)];

cond=[1;0;0;0;0;0];

options=odeset('RelTol',1e-8,'OutputFcn',@odewbar);
[t,y]=ode15s(dxdt,[0,1e-7],cond,options);

figure(1)
plot(t,y(:,1),'k',t,y(:,4),'b',t,y(:,6),'r')
legend('$\rho_{gg}$','$\rho_{kk}$','$\rho_{ee}$','interpreter','latex')
set(gca,'GridAlpha',0.7,'MinorGridAlpha',0.7,'FontSize',28,...
    'FontName','Times')
set(findobj(gcf,'Type','Line'),'LineWidth',1.5)
grid on
grid minor