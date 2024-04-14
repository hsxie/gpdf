% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-05-28
% test zetaph.m 2d, using zetaph.m

close all; clear; clc;
a=2;
[xx,yy]=meshgrid(-a:0.05*a:a,-a:0.05*a:a);
z=xx+1i*yy;
[Zp2,Z2]=zetaph(z);

h=figure('unit','normalized','Position',[0.01 0.37 0.6 0.55]);
set(gcf,'DefaultAxesFontSize',15);

subplot(221);surf(xx,yy,real(Zp2)); xlabel('x'); ylabel('y'); 
title('(a) Re(Z''), surface plot'); axis tight;
subplot(222);surf(xx,yy,imag(Zp2)); xlabel('x'); ylabel('y'); 
title('(b) Im(Z''), surface plot'); axis tight;

a=6; 
x=-a:0.01*a:a;y=-a:0.01*a:a;
[xx,yy]=meshgrid(x,y);
z=xx+1i*yy;
[Zp2,Z2]=zetaph(z);

subplot(223);
imagesc(x,y,real(Z2));axis xy square; caxis([-1 1]); colorbar;
xlabel('x'); ylabel('y'); title('(c) Re(Z), scales image contour plot');
% subplot(223);pcolor(xx,yy,real(Z2)); shading interp;
subplot(224);
contourf(xx,yy,imag(Z2));axis xy square; colorbar;
xlabel('x'); ylabel('y'); title('(d) Im(Z), filled contour plot');
% subplot(224);pcolor(xx,yy,imag(Z2)); shading interp;


