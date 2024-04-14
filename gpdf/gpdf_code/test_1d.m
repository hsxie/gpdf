% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-05-28
% test zetaph.m 1d, using zetaph.m

close all; clear; clc;

z0=-10:0.1:10;
z=z0-0.00i;
% z0=(-1:0.01:1);
% z=z0.*1i-1.00;

% faddeeva.m in real line is calculate using erfc()
Z1=faddeeva(z)*1i*sqrt(pi); 
Zp1=-2*(1+z.*Z1);
[Zp2,Z2]=zetaph(z);

h=figure('unit','normalized','Position',[0.01 0.47 0.6 0.45]);
set(gcf,'DefaultAxesFontSize',15);


subplot(121);
plot(z0,real(Z1),'g--',z0,real(Z2),'r.',z0,real(Zp1),...
    'm--',z0,real(Zp2),'b.','Linewidth',2);
xlabel('Re(z)'); legend('exact Re(Z)','GPDF Re(Z)',...
    'exact Re(Z'')','GPDF Re(Z'')');
legend('boxoff'); ylim([-2.0,2]);
title(['(a) Re(Z) and Re(Z''), Im(z)=',num2str(imag(z(1)))]);

subplot(122);
plot(z0,imag(Z1),'g--',z0,imag(Z2),'r.',z0,imag(Zp1),...
    'm--',z0,imag(Zp2),'b.','Linewidth',2);
xlabel('Re(z)');legend('exact Im(Z)','GPDF Im(Z)',...
    'exact Im(Z'')','GPDF Im(Z'')');
legend('boxoff'); ylim([-1.6,2.5]);
title(['(b) Im(Z) and Im(Z''), Im(z)=',num2str(imag(z(1)))]);

