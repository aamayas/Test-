clear;
%%ID script to draw figure 8
%% for no rotation and simple basin, eta can be found analytically  
function [u,v,w]=uvw(eta_x,eta_y,z,h,h_y,xmax);
  nz=length(z);
  u=zeros(size(z));v=u;w=u;
  for mz=1:nz
    if(z(mz)>=-h)
      u(mz)=eta_x*(z(mz)^2-h^2)/2+z(mz)+h;%   eqn 1.15
      v(mz)=eta_y*(z(mz)^2-h^2)/2;%    eqn 1.15
      w(mz)=h_y*eta_y*z(mz)*(z(mz)^2-h^2)/2/h;%   equation 1.33
    end
  end
endfunction

%%  MAIN
%%set up plots
set(gcf,'position',[100,100,1000,400]);clf

%% compute eta from analytical solution
tau=1;
xmax=3;
nx=101;
ny=18;
nz=31;
x=xmax*linspace(-1,1,nx);dx=diff(x(1:2));
y=linspace(-1,1,ny);dy=diff(y(1:2));
z=linspace(0,-1,nz);
h=0.01+0.99*(1-abs(y));
h_y=firstder(h,dy);
%% for no rotation eta can be found analytically: 
eta=etaana(tau,x,h);
[eta_x,eta_y]=mygradient(eta,dx,dy);
%% drawing two sections close to the closed end at locations ix
for isec=1:2
  u=v=w=zeros(ny,nz);
  ix=92+8*(isec-1);
  %% have chosen ny points along y-axis
  for ky=1:ny
    %% find u,v,w profiles at eah y location
    [u(ky,:),v(ky,:),w(ky,:)]=uvw(eta_x(ky,ix),eta_y(ky,ix),z,h(ky),h_y(ky),xmax);
  end

  axes('position',[0.09+0.45*(isec-1),0.125,0.43,0.85]);
  colormap([.75 .75 .75;1 1 1]);% gray=.75, white=1
  contourf(y,z,u',2.5e-2*[-10:2:10]);hold on

  set(gca,'xdir','reverse');% the model y axis is ploted along the plot x axis
  xlabel('y');
  if isec==1;
    ylabel('z');% label y axis on left plot only
  else
    set(gca,'ytickmode','manual','ytick',[],'yticklabelmode','manual','yticklabel',[]);
  endif
  
  %%draw vectors for v,w
  for ky=1:ny
    for kz=1:2:nz
      if z(kz)>-h(ky)
        Arrow(y(ky),z(kz),v(ky,kz),w(ky,kz),1.8,2,'k',0.075);
      endif
    endfor
  endfor
  fill([1.,1.,0.],[0.,-1.,-1],'w');
  fill([-1.,-1.,0.],[0.,-1.,-1],'w');
endfor
set(gcf,'paperpositionmode','manual','paperposition',[0.25 0.25 10 4])
%print -dpng ../../words/figs/C1/figure6.png;% uncomment to save
