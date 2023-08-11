clear;
%%ID show psi in IdealBasin with delta=1/6  script to draw figure 3
%%ID
load XTHNoTaper.mat
delta=1/6;%input('delta');
theta=0;%input('theta radians');% wind parallel to x-axis
%% need function psinum.m to be in this directory
psi=psinum(theta,delta,xp,tp,h,Shore);
%% save output for use in drawing local velocity sections.
save psinum.mat theta delta xp tp h psi
set(gcf,'position',[900,800,1200,450]);clf;
set(0, 'DefaultAxesFontsize',16)
set(0, 'DefaultTextFontsize',16)

axes('position',[0.08,0.1,0.85,0.85]);
tricontour(tp,xp(:,1),xp(:,2),psi,2.5e-3*[-10:10],'k',1);
%% everything below here is to make the plot look good
hold on
axis([min(xp(:,1)),max(xp(:,1)),min(xp(:,2)),max(xp(:,2))+.3]);
axis equal
plot([-3,3,3,-3,-3],[-1,-1,1,1,-1],'k','linewidth',2);% basin frame

%% wind stress arrow:
Arrow(-0.5,1.2,1,0,1,2,'k',0.25);
text(-.75,1.2,'\tau','fontsize',22,'fontweight','bold')
%% circulation arrows:
Arrow(.25,0,-.5,0,1,2,'k',0.2);% vector along x-axis toward -x
Arrow(-.25,-0.71,.5,0,1,2,'k',0.2);% vector  toward -x
Arrow(-.25,0.71,.5,0,1,2,'k',0.2);% vector  toward x
%% locate sections
plot([0,0],[-1,1],'k--','linewidth',1);
plot([2.9,2.9],[-1,1],'k-.','linewidth',1);
xlabel('x','fontsize',18,'fontweight','bold')
ylabel('y','fontsize',18,'fontweight','bold')
box off
set(gcf,'paperpositionmode','manual','paperposition',[0.25 0.25 12 4.5])
%print -dpng ../../words/figs/C2/figure3.png
