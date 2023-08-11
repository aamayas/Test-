clear;

function[x,y]=latlon2xy(lon,lat,minlon,minlat);
  %% minlon and minlat determine origin of x,y coordinates
  %% 1 degrre of latitude ~ 111 km
  theta=min(lat)*pi/180;
  x=111*cos(theta)*(lon-minlon);
  y=111*(lat-minlat);
endfunction
  
%%ID 2 figures: lhs BC bathy; rhs psi
%%
set(gcf,'position',[200,1200,1100,600]);
clf;
%% Left: plot bathymetry, available in XTHBC.mat
load XTHBC.mat
axes('position',[0.1,0.13,0.35,0.85])
plot(xcoast,ycoast,'k','linewidth',2);
hold on;
tricontour(tp,xp(:,1),xp(:,2),h,[0:5:35],'k',1);
axis([0,35,0,40])
axis equal;
%% the following 12 lines to label the x axis in longitudes and latitudes
LONTICK=[-111.95:.1:-111.65];% want ticks at these lon;
XTICK=111*cos(26.5*pi/180)*(LONTICK'+111.97);
LONTICK=num2str(transpose(LONTICK));
LONcell=cellstr(LONTICK);
set(gca,'xtickmode','manual','xtick',XTICK,'xticklabelmode','manual','xticklabel',LONcell);
xlabel('longitude','fontsize',14','fontweight','bold');
LATTICK=[26.5:.1:26.9];% want ticks at these lon;
YTICK=111*(LATTICK'-26.539);
LATTICK=num2str(transpose(LATTICK));
LATcell=cellstr(LATTICK);
set(gca,'ytickmode','manual','ytick',YTICK,'yticklabelmode','manual','yticklabel',LATcell);
ylabel('latitude','fontsize',14','fontweight','bold');
title('Coastline and Bathymtry','fontsize',18,'fontweight','bold');
drawnow

%% Right: plot transport streamlines on RHS
axes('position',[0.6,0.13,0.35,0.85])
%% begin by redrawing coast, but this time using x,y axes (see BuildXTH)
plot(xcoast,ycoast,'k','linewidth',2);
hold on
h=h/max(h);%h is now relative depth (as requd by psinum)
theta=-40*pi/180.;%% wind orientation;


%% change 9 June 2012: was
%%psi=psinumf0(theta,xp,tp,h);
delta=100;
psi=psinum(theta,delta,xp,tp,h,Shore);
%% represent and label tau
%%Arrow2d([25.,25.;25.+10*cos(theta),25.+10*sin(theta)],'k',0.25);
Arrow(25.,25.,10*cos(theta),10*sin(theta),1,1.5,'k',2);
text(26,22,'\tau','fontsize',18,'fontweight','bold');% label wind stress
tricontour(tp,xp(:,1),xp(:,2),psi,[0 0],'k',2);
tricontour(tp,xp(:,1),xp(:,2),psi,5.e-2*[-10:10],'k',1);
%% draw direction of psi circulation

Arrow(20.818,15.909,-1,2,2,1.5,'k',1.5);% draw return current
Arrow(19.909,13.909,1,-2,1,1.5,'k',1);
Arrow(24.818,16.455,2,-2,2,1.5,'k',1);
Arrow(14.273,35.545,.2,-2,2,1.5,'k',1);

axis([0,35,0,40])
axis equal
%plot(xp(basinedge,1),xp(basinedge,2),'k*');% plot coastline
xlabel('x (km)','fontsize',14','fontweight','bold');
ylabel('y (km)','fontsize',14','fontweight','bold');
title('Tranport Streamfunction \psi','fontsize',18,'fontweight','bold');

set(gcf,'paperpositionmode','manual','paperposition',[0.25 0.25 11 6])
%%print -dpng ../../words/figs/C1/figure7.png
