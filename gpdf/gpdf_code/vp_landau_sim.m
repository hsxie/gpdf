% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-11-16 15:50
% half spectral method ([Xie2012] arXiv:1211.5984) for Landau damping 
% simulation
% ion immobile, e=1, m=1, ...
% d(df(v,t))/dt=-ikv*(df)+(dE)*d(f0)/dv
% ik(dE(t))=-\int{df}dv
% 2012-11-27: exp(-v^2) to exp(-v^2/2), Landau damping match theory now!
% 2013-02-22: 1/t damping Case-Van Kampen mode ([Xie2013a] arXiv:1304.5883)
% 2013-05-28: used for benchmark GPDF (Generalized Plasma Dispersion
%             Function), [Xie2013b] arXiv:1305.6476
close all; clear; clc;

% Landau damping data
data = [0.1   1.0152   -4.75613E-15
0.2   1.06398   -5.51074E-05
0.3   1.15985   -0.0126204
0.4   1.28506   -0.066128
0.5   1.41566   -0.153359
0.6   1.54571   -0.26411
0.7   1.67387   -0.392401
0.8   1.7999   -0.534552
0.9   1.92387   -0.688109
1.0   2.0459   -0.85133
1.5   2.63233   -1.77571
2   3.18914   -2.8272];
id=2;
k=data(id,1); wr=data(id,2); wi=data(id,3);
k=0.15;

nv=256*8; dt=0.01; nt=2000;
vmax=20; vmin=-vmax;
dv=(vmax-vmin)/nv; TR=2*pi/(k*dv);
vv=vmin:dv:vmax;

% % initial
% df=0.*vv+0.1.*exp(-(vv-1.0).^2); % for some distributions, sensitive here
df=0.*vv+0.1.*exp(-(vv-0.0).^2);
dE=0.01;
% % dE=(1i/k)*sum(df)*dv;

Fn=2; % choose initial distribution function

if(Fn==1)
    % Maxwellian
    vd=0.0;
    f0=exp(-(vv-vd).^2./2)/(sqrt(2*pi)); % f0=exp(-vv.^2), Maxwellian
    df0dv=-(vv-vd).*exp(-(vv-vd).^2./2)/(sqrt(2*pi));
%     % incomplete Maxwellian
%     nu=-0.1*sqrt(2);
%     if(nu>-vmax && nu<vmax)
%         ind=find(vv<=nu);
%         df0dv(ind)=0;
%     %     df0dv(ind(end))=(f0(ind(end)))/dv;
%     end
elseif(Fn==2)
    % Kappa    
    kappa=1;
    T=1; vt=sqrt(2*T);
    if(kappa==1)
        wr=1; wi=-k*vt;
    end
    Akappa=(gamma(kappa)/gamma(kappa-0.5)/(sqrt(pi*kappa)*vt));
    % f0=Akappa.*(1+vv.^2/(kappa*vt^2)).^(-kappa); % Kappa
    df0dv=-Akappa.*(2.*vv)./(vt.^2.*(vv.^2./(kappa.*vt.^2)+1).^(kappa+1));
elseif(Fn==3)
    % Triang Fp='((v>-1)&(v<=0))-((v>0)&(v<=1))';
    df0dv=((vv>-1)&(vv<0))-((vv>0)&(vv<1));
else
    % Slowing down
    vt=sqrt(2);
    df0dv=(((vv>=0)&(vv<=4)).*(-3.*vv.^2)./(vv.^3+vt^3).^2+...
        ((vv<0)&(vv>-4)).*(3.*vv.^2)./((-vv).^3+...
        vt^3).^2)*3*sqrt(3)*vt^2/(4*pi);
end

dEt=zeros(1,nt+1);
dEt(1)=dE;
tt=linspace(0,nt*dt,nt+1);
 
h=figure('unit','normalized','Position',[0.01 0.57 0.6 0.35]);
set(gcf,'DefaultAxesFontSize',15);

% plot(vv,df0dv,'r--');

for it=1:nt
    % RK-4, 1st step
    pdf1=-1i*k.*vv.*df+dE.*df0dv;
    % RK-4, 2nd step
    dftmp=df+0.5.*dt.*pdf1;
    dEtmp=(1i/k)*sum(dftmp)*dv;
    pdf2=-1i*k.*vv.*dftmp+dEtmp.*df0dv;
    % RK-4, 3rd step
    dftmp=df+0.5.*dt.*pdf2;
    dEtmp=(1i/k)*sum(dftmp)*dv;
    pdf3=-1i*k.*vv.*dftmp+dEtmp.*df0dv;
    % RK-4, 4th step
    dftmp=df+dt.*pdf3;
    dEtmp=(1i/k)*sum(dftmp)*dv;
    pdf4=-1i*k.*vv.*dftmp+dEtmp.*df0dv;
    % RK-4, push
    df=df+dt./6.0.*(pdf1+2.0.*pdf2+2.0.*pdf3+pdf4);
    dE=(1i/k)*sum(df)*dv;
    
    dEt(it+1)=dE;
    
    % plot
    if(mod(it-1,nt/20)==0 || it==nt)
        subplot(121);
        plot(tt(2:it),real(dEt(2:it)),tt(2:it),...
            imag(dEt(2:it)),'-.','LineWidth',2);
        title(['(a) k=',num2str(k),', dv=',num2str(dv),...
            ', dt=',num2str(dt)]); grid on;
        legend('Re[\delta E]','Im[\delta E]',2); legend('boxoff');
        xlabel('t'); ylabel('\delta E'); xlim([0,tt(end)]);
        ttmp=tt(1:it);lndE=log(real((dEt(1:it))));
        lndEe=log(abs(real(dEt(2).*exp((-1i*wr+wi).*ttmp))))-15;
        if(dt*it<=100)
            subplot(122);title(['(b) \omega_{theory}=',num2str(wr+1i*wi)]);
            plot(tt(1:it),lndEe,'r:','LineWidth',2);
        end
        subplot(122);plot(tt(1:it),real(lndE),'LineWidth',2);hold on;
        xlabel('t');ylabel('ln|\delta E|'); xlim([0,tt(end)]);
        legend('simulation','theory',3);legend('boxoff');
       
        pause(0.01);
    end
end

% Find the corresponding indexes of the extreme max values 
it0=floor(nt*1/20); it1=floor(nt*4/5);
yy=lndE(it0:it1);
extrMaxIndex = find(diff(sign(diff(yy)))==-2)+1;
subplot(122);hold on;
t1=tt(it0+extrMaxIndex(1));t2=tt(it0+extrMaxIndex(end));
y1=yy(extrMaxIndex(1));y2=yy(extrMaxIndex(end));
plot([t1,t2],[y1,y2],'r*--','LineWidth',2);
omega=pi/((t2-t1)/(length(extrMaxIndex)-1));
gammas=(real(y2)-real(y1))/(t2-t1);
title(['(b) \omega^S=',num2str(omega),', \gamma^S=',num2str(gammas)]);

% str=['vmax=',num2str(vmax),',nv=',num2str(nv),...
%     ',dt=',num2str(dt),',k=',num2str(k)];
% print('-dpng',[str,'.png']);
% print('-dpsc',[str,'.eps']);
