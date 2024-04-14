% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-05-28
% root finding for Langmuir wave, using zetaph.m

close all; clear; clc;
zetap=@(x)zetaph(x);
f=@(x,k)k*k-zetap(x)/2;
% f=@(x,k)k*k-zetap(x);
w=[];
% k=[0.02:0.02:0.2,0.2:0.1:1.0,1.5:0.5:10];
k=[0.1:0.1:1.0,1.5:0.5:5];
xg=1.5e-0-1i*1e-3;
for kk=k
    options=optimset('Display','off');
    wg=xg/(sqrt(2)*kk);
    x=fsolve(f,wg,options,kk)*sqrt(2)*kk;
%     wg=xg/kk;
%     x=fsolve(f,wg,options,kk)*kk;
    w=[w,x];
    xg=x;
end
[k,ind]=sort(k);
w=w(ind);
wre=real(w); wie=imag(w);
wrt=1.0+1.5.*k.*k;
wit=-sqrt(pi/8).*exp(-1.0./(2.0.*k.^2)-1.5)./(k.^3);

h=figure('unit','normalized','Position',[0.01 0.57 0.65 0.35]);
set(gcf,'DefaultAxesFontSize',15);

subplot(121);
% loglog(k,wre,'r',k,wrt,'b--','LineWidth',2);
% loglog(k,abs(wre),'r','LineWidth',2); hold on;
semilogx(k,abs(wre),'r+','LineWidth',2); hold on;
xlabel('k\lambda_D');ylabel('\omega/\omega_p');
xlim([min(k),max(k)]); 
ylim([0.0,1e1]);
title('(a) real frequency v.s. k');

subplot(122);
% loglog(k,-wie,'+r--',k,-wit,'b--','LineWidth',2);
% semilogx(k,wie,'r',k,wit,'b--','LineWidth',2);
semilogx(k,-wie,'r+','LineWidth',2); hold on;
xlabel('k\lambda_D');
ylabel('-\gamma/\omega_p');
xlim([min(k),max(k)]);
% ylim([-10.0,2]);
title('(b) damping rate v.s. k');




