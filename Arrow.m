
%% %%ID usage: Arrow(x0.y0,u,v,scale,width,color,head_size)
%% x0,y0 =arrow origin
%% u,v,vector components
%% width linewidth for shaft
%% head_size size of tip
%% test program
%%color=[0,.5,0];scale=1;width=1;head_size=.1;
%%x0=y0=0;
%%axis equal
%%for kt=1:16
%%  v=sin(2*pi*(kt-1)/15);
%%  u=cos(2*pi*(kt-1)/15);
%%  Arrow(x0,y0,u,v,scale,width,color,head_size)
%%  hold on
%%endfor

function Arrow(x0,y0,u,v,scale,width,color,head_size);
  %% at this stage I need axis equal, otherwise distorted
  %% to get AR:
  xl = get(gca,'xlim');
  yl = get(gca,'ylim');
  dx = xl(2)-xl(1);
  dy = yl(2)-yl(1);
  AR = dy/dx ;
  %% draw the shaft
  u0 = scale*u ;
  v0 = scale*v ;
  xx = x0 + [0 u0] ;
  yy = y0 + [0 v0] ;
  % draw shaft
  plot(xx,yy,'color',color,'linewidth',width)
  %% add the tip xa,ya: coordinates of unrotated tip-triagle apex
  %% values 1 and 0.3 set the arrow shape
  xa=[-1,0,-1]*head_size;ya=[-0.3,0,0.3]*head_size;
  theta=atan2(v,u);%atan2(y2-y1,x2-x1);
  R=[cos(theta),-sin(theta);sin(theta),cos(theta)];
  %% rotate and translate in next line
  rot=R*[xa;ya];
  %% next line makes tip size proportional to vector length??
  xr=x0+u0+rot(1,:);yr=y0+v0+rot(2,:);
  patch(xr,yr,'edgecolor',color,'facecolor',color)
endfunction

