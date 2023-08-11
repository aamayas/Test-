clear;
%%ID Show effect of wind angle on psi for delta=1/6

delta=1/6;

%% load XTH for Bahia Concepcion
load XTHBC.mat
h=h/max(h);

%% get two solutions: psix wind along x; psiy wind along y
%%this way I can construct solutions for any theta
psix=psinum(0,delta,xp,tp,h,Shore);
psiy=psinum(pi/2,delta,xp,tp,h,Shore);%% theta in radians

%% plot 10 solutions

set(gcf,'position',[10,1000,1000,500]);clf;
for klr=1:5
  thetalr=pi/180*10*(klr-1);
  for kud=1:2
    theta=thetalr+pi/180*50*(kud-1);
    theta=-theta;
    psi=psix*cos(theta)+psiy*sin(theta);
    axes('position',[-.01+(klr-1)/5,0.52-(kud-1)/2,0.25,0.48]);
    tricontour(tp,xp(:,1),xp(:,2),psi,[0 0],'k',2);hold on
    tricontour(tp,xp(:,1),xp(:,2),psi,1e-2*[-12:12],'k',1);
    Arrow(25.,25.,cos(theta),sin(theta),5,1.5,'k',1.5);% draw wind stress
    text(22.,25,'\tau','fontsize',20,'fontweight','bold');% label wind stress
    axis([min(xp(:,1)),max(xp(:,1)),min(xp(:,2)),max(xp(:,2))]);
    axis equal
    set(gca,'xtickmode','manual','xtick',[],'xticklabelmode','manual','xticklabel',[]);
    set(gca,'ytickmode','manual','ytick',[],'yticklabelmode','manual','yticklabel',[]);
    axis off;
    box off;
    drawnow;
  endfor
endfor

%% to print
set(gcf,'paperpositionmode','manual','paperposition',[0.25 0.25 10 5])
%print -dpng ../../words/figs/C2/f6PsiofthetaBC.png
