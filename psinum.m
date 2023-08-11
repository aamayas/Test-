function psi=psinum(theta,delta,xp,tp,h,Shore);
  %% solves equation 1.20 for a general domain using FEM
  %% for delta < infnty

  tausupx=cos(theta);
  tausupy=sin(theta);
  h=h/max(h);
  
  np=length(xp(:,1));nt=length(tp(:,1));
  Interior=setdiff([1:np]',Shore);
  if delta > 10;% this is a case where rotation can be ignored
    A=1./h.^3;
    B=0;
    C=1./h.^2;
    D=0;
  else
    %% compute coefficients A, B, C, D in PDE
    alpha=(1+i)/delta;
    alphah=alpha*h;
    cah=cosh(alphah);tah=tanh(alphah);
    tau=(cah-1)./cah/alpha^2;Qtau=real(tau);Rtau=imag(tau);
    eta=(tah-alphah)./alpha^3;Qeta=real(eta);Reta=imag(eta);
    den=Qeta.^2+Reta.^2;
    A=-Qeta./den;
    B=-Reta./den;
    C=-(Qeta.*Qtau+Reta.*Rtau)./den;
    %%  CAREFUL DO NOT CONFUSE THIS D WITH GRADIENT OPERATORS Dx and Dy  %%
    D=-(Qeta.*Rtau-Reta.*Qtau)./den;
  endif
  %% Initialize global matrices
  M=Dx=Dy=PDE=sparse(np,np);
  load QuadW12; % load Quadrature points and weights
  for kt=1:nt;% element by element construction og global matrices
    nodes=tp(kt,:);ind=tp(kt,:);% they are the same for linera test functions
    %% next lines finds element matrices Ke, Me, Dex, Dey
    %% and forcing column vector Fe 
    [Cs,a,b,c,theta]=LinearElement(xp(nodes,:));
    Me = ElementMatrixMe (Cs,a,b,c);
    [Dex,Dey] = ElementMatrixDexy (Cs,a,b,c,theta);
    %% PDE has variable coefficient A(h): evaluate Ke by quadrature
    x1=xp(nodes(1),1);x2=xp(nodes(2),1);x3=xp(nodes(3),1);
    y1=xp(nodes(1),2);y2=xp(nodes(2),2);y3=xp(nodes(3),2);
    detF=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
    Finv=[y3-y1,x1-x3;y1-y2,x2-x1]/detF;
    Area=detF/2;
    phi_X=[-1,-1;1,0;0,1];% test function derivative in canonical traingle
    phi_x=phi_X*Finv;% test function derivative here
    Ke=zeros(3,3);% Initialize Element Matrix Ke
    for kq=1:length(wt); % 1st quadrature loop
      %% next line: % A evaluated at xq,yq
      Av=A(nodes(1))*(1-xq(kq)-yq(kq))+A(nodes(2))*xq(kq)+A(nodes(3))*yq(kq);
      if delta <=10
	Bv=B(nodes(1))*(1-xq(kq)-yq(kq))+B(nodes(2))*xq(kq)+B(nodes(3))*yq(kq);
      endif
      for m=1:3
	for n=1:3
	  if delta <=10
	  Ke(m,n)=Ke(m,n)+Area*wt(kq)*phi_x(m,1)*(Av*phi_x(n,1)+Bv*phi_x(n,2));
	  Ke(m,n)=Ke(m,n)-Area*wt(kq)*phi_x(m,2)*(Bv*phi_x(n,1)-Av*phi_x(n,2));
	  else
	    Ke(m,n)=Ke(m,n)+Area*wt(kq)*phi_x(m,1)*Av*phi_x(n,1);
	    Ke(m,n)=Ke(m,n)+Area*wt(kq)*phi_x(m,2)*Av*phi_x(n,2);
	  endif	    
	endfor
      endfor
    endfor; % end of 1st quadrature loop
    %% next lines assemble global matrices from element matrices.
    M(nodes,nodes)=M(nodes,nodes)+Me;
    Dx(nodes,nodes)=Dx(nodes,nodes)+Dex;
    Dy(nodes,nodes)=Dy(nodes,nodes)+Dey;
    PDE(ind,ind)=PDE(ind,ind)-Ke;
  endfor
  %% had to wait for M, Dx and Dy to be assembled for this next step
  if delta <=10
    RHS=M\(tausupx*(Dx*D-Dy*C)+tausupy*(Dx*C+Dy*D));
  else
    RHS=M\(tausupy*Dx*C-tausupx*Dy*C);
  endif
  F=zeros(np,1);
  for kt=1:nt;% element by element construction og global matrices
    nodes=tp(kt,:);ind=tp(kt,:);
    %% element vector Fe and global vector F assembled separately
    %% also using quadrature
    x1=xp(nodes(1),1);x2=xp(nodes(2),1);x3=xp(nodes(3),1);
    y1=xp(nodes(1),2);y2=xp(nodes(2),2);y3=xp(nodes(3),2);
    detF=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
    Finv=[y3-y1,x1-x3;y1-y2,x2-x1]/detF;
    Area=detF/2;
    Fe=zeros(3,1);% initialize forcing element vector
    for kq=1:length(wt); %2nd quadrature loop
      %% value of RHS at this quadrature point
      RHSv=RHS(nodes(1))*(1-xq(kq)-yq(kq))+RHS(nodes(2))*xq(kq) ...
	   +RHS(nodes(3))*yq(kq);
      phi(1)=1-xq(kq)-yq(kq);% value of test function 1 at this point
      phi(2)=xq(kq);
      phi(3)=yq(kq);
      for m=1:3
	Fe(m)=Fe(m)+Area*wt(kq)*RHSv*phi(m);        
      endfor
    endfor; % end of2nd  quadrature loop
    F(ind)=F(ind)+Fe;% assemble in to global forcing vector F
  endfor %% end of 2nd assembly loop

  %%set boundary conditions on psi=0 on periphery of basin
  psi=zeros(np,1);% intialize all psi to zero (including boundary points)

  %% solve the linear system, values on boundaries already known
  psi(Interior)=PDE(Interior,Interior)\F(Interior);

endfunction


