clear;
%%ID Draw a section across a basin, show velocity normal to the section
%%ID with contours, and velocities in the section as arrows.
%%ID This is for the solution currently in saved file psinum.mat
%%ID

load psinum.mat

nz=16;%input('no of levels on deepest profile');
%% in the next line I choose z to decrease from the surface,
z=linspace(0,-1,nz);

np=length(xp(:,1));ne=length(tp(:,1));

%% Build FEM global matrices for Mass (m) X-grad (Dx) Y-grad (Dy)
M=Dx=Dy=sparse(np,np);% Initialize Global Matrices as SPARSE
for ke=1:ne
  nodes=tp(ke,:);
  [Cs,a,b,c,alpha]=LinearElement(xp(nodes,:));
  Me = ElementMatrixMe (Cs,a,b,c);
  M(nodes,nodes)=M(nodes,nodes)+Me;
  [Dex,Dey] = ElementMatrixDexy (Cs,a,b,c,alpha);
  Dx(nodes,nodes)=Dx(nodes,nodes)+Dex;
  Dy(nodes,nodes)=Dy(nodes,nodes)+Dey;
endfor
%% find eta_x and eta_yc
psi_x=M\(Dx*psi);psi_y=M\(Dy*psi);
%% compute coefficients A, B, C, D in PDE
gamma=(1+i)/delta;
gammah=gamma*h;
cah=cosh(gammah);tah=tanh(gammah);
tau=(cah-1)./cah/gamma^2;Qtau=real(tau);Rtau=imag(tau);
eta=(tah-gammah)./gamma^3;Qeta=real(eta);Reta=imag(eta);
den=Qeta.^2+Reta.^2;
A=-Qeta./den;
B=-Reta./den;
C=-(Qeta.*Qtau+Reta.*Rtau)./den;
%%  CAREFUL DO NOT CONFUSE THIS D WITH GRADIENT OPERATORS Dx and Dy  %%
D=-(Qeta.*Rtau-Reta.*Qtau)./den;
%% wind stress as compnents, vector length-1
taux=cos(theta);tauy=sin(theta);
%% apply equation 1.21
eta_x=A.*psi_y-B.*psi_x+C*taux-D*tauy;
eta_y=-A.*psi_x-B.*psi_y+C*tauy+D*taux;
%% comput u,v,w at each node

gamma=(1+i)/delta;% always the same
for kp=1:np
  gammah=gamma*h(kp);%  depth is negative
  cah=cosh(gammah);
    for kz=1:nz; % each depth
    if z(kz) >=-h(kp)
      qirtau=sinh(gamma*(z(kz)+h(kp)))/cah/gamma;% (2.7)
      qireta=(1-cosh(gamma*z(kz))/cah)/gamma^2;%eqn 1.13
      u(kp,kz)=real(qirtau*(taux+i*tauy)-(eta_x(kp)+i*eta_y(kp))*qireta);% (2.5)
      v(kp,kz)=imag(qirtau*(taux+i*tauy)-(eta_x(kp)+i*eta_y(kp))*qireta);% (2.6)
    else
      u(kp,kz)=0;
      v(kp,kz)=0;
    endif
  endfor
endfor

%% w: make v_y(kp,kz)
for kz=1:nz
  u_x(:,kz)=M\(Dx*u(:,kz));
  v_y(:,kz)=M\(Dy*v(:,kz));
endfor
dz=1/(nz-1);
for kp=1:np
  %% 2 minus signs in next line to remind that we are integrating\
  %% -du/dx-dv/dy from 0 (where w=0) to -1 (so that dz<0). 
  w(kp,:)=integral(-(u_x(kp,:)+v_y(kp,:)),-dz);
endfor

%% MIDPOINT: Local velocities u, v, w are available at all xp, for z>-h

%% SPECIAL for figure4.m: Draw 2 sections side by side
set(gcf,'position',[100,100,1800,600]);clf
colormap([.75 .75 .75;1 1 1]);% gray=.75, white=1
%% Draw contour and arrows for a section located at x

%% in this special grid, at some value of 1<kx<=nx
% the list of positions will be kx:nx:np-(nx-kx)
%% we are looking along + x, meaning that the
%% y-axis needs to be reversed (-1 on right) 
%% masking off areas beneath the bottom with fill
ny=61;%input('how many profiles in a section?');

for ksec=1:2
  x=0;% for ksec==1
  if ksec == 2; x=2.9; endif
  %% interpolate from xp points to section
  x_S=x*ones(1,ny);
  y_S=linspace(-1,1,ny);
  h_S=griddata(xp(:,1),xp(:,2),h,x_S,y_S);
  u_S=v_S=w_S=zeros(ny,nz);
  for kz=1:nz
    u_S(:,kz)=griddata(xp(:,1),xp(:,2),u(:,kz),x_S,y_S);
    v_S(:,kz)=griddata(xp(:,1),xp(:,2),v(:,kz),x_S,y_S);
    w_S(:,kz)=griddata(xp(:,1),xp(:,2),w(:,kz),x_S,y_S);
  endfor
  
  %% show u as contour and v,w as vectors
  axes('position',[0.08+0.424*(ksec-1),0.1,0.4,0.85]);
  contourf(y_S,z,u_S',2.5e-2*[-10:2:10]);hold on
  caxis([-2*eps,-eps]);% orces transition between gray and white to be 0
  set(gca,'xdir','reverse');% the model y axis is ploted along the plot x axis
  
  xlabel('y','fontsize',18,'fontweight','bold');
  if ksec == 1;ylabel('z','fontsize',18,'fontweight','bold');endif;
  if ksec == 2;set(gca,'ytickmode','manual','ytick',[],'yticklabelmode','manual','yticklabel',[]);endif;
  %%draw vectors for v,w
  for ky=1:3:ny
    for kz=1:nz
      if z(kz)>-h_S(ky)
        Arrow(y_S(ky),z(kz),v_S(ky,kz),w_S(ky,kz),1.8,2,'k',0.05);
      endif
    endfor
  endfor
  fill([1.,1.,0.],[0.,-1.,-1],'w');
  fill([-1.,-1.,0.],[0.,-1.,-1],'w');
endfor
set(gcf,'paperpositionmode','manual','paperposition',[0.25 0.25 18 6])
%%print -dpng ../../words/figs/C2/figure4.png

