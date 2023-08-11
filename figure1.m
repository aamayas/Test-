%%ID compute and plot q(z) and r(z) from equation 2.7
%%ID draws figure 2.1

clear
delta=0.2;
nz=51;
z=linspace(-1,0,nz);

alpha=(1+i)/delta;
qirtau=sinh(alpha*(z+1))/alpha/cosh(alpha);
qireta=(cosh(alpha*z)/cosh(alpha)-1)/alpha^2;

set(gcf,'position',[500,1000,1000,700]);clf
axes('position',[0.1,0.1,0.35,0.85]);
plot(real(qirtau),z,'linewidth',2);
hold on
plot(imag(qirtau),z,'linewidth',2);
xlabel('q^\tau, r^\tau');
legend('q^\tau','r^\tau','location','southeast')
axes('position',[0.55,0.1,0.35,0.85]);
plot(real(qireta),z,'linewidth',2);
hold on
plot(imag(qireta),z,'linewidth',2);
legend('q^\eta','r^\eta','location','southeast')
xlabel('q^\eta, r^\eta');
ylabel('z')
%%print -dpng ../../words/figs/C2/figure1.png
