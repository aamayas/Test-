clear;
%%ID solution for psi(y) and [u]=-d psi/dy
function [A,B,C,D]=ABCD(delta,h);
  %% evaluate functions A,B.C.d feined in equation 2.15
  alpha=(1+i)/delta;
  alphah=alpha*h;
  cah=cosh(alphah);tah=tanh(alphah);
  tau=(cah-1)/cah/alpha^2;Qtau=real(tau);Rtau=imag(tau);
  eta=(tah-alphah)/alpha^3;Qeta=real(eta);Reta=imag(eta);
  den=Qeta^2+Reta^2;
  A=-Qeta/den;
  B=-Reta/den;
  C=-(Qeta*Qtau+Reta*Rtau)/den;
  D=-(Qeta*Rtau-Reta*Qtau)/den;
endfunction

%% MAIN Section 2.4    makes figure 2.2

clf;
ny=51;y=linspace(-1,1,ny);h=0.05+0.95*(1-abs(y));dy=y(2)-y(1);

delta=1/6; %% was 0.1
sumoneoA=0;sumCoA=0;
for ky=1:ny
[A(ky),B(ky),C(ky),D(ky)]=ABCD(delta,h(ky));
sumoneoA=sumoneoA+1/A(ky);
sumCoA=sumCoA+C(ky)/A(ky);
end
for ky=1:ny
psi_y(ky)=(sumCoA/sumoneoA-C(ky))/A(ky);%% equation 2.20
end

plot(y,-psi_y,'k--','linewidth',2)
hold on;
xlabel('y','fontsize',16')
ylabel('\psi, [u]','fontsize',16);
psi=integral(psi_y,dy);
plot(y,psi,'k','linewidth',2);
plot([-1 1],[0 0],'k')
set(gca,'xdir','reverse');% the model y axis is ploted along the plot x axis
%print -dpng ../../words/figs/C2/figure2.png
