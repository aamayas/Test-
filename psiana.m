function psi=psiana(tau,x,h);
  %%ID script to compute Mathieu 2002 solution for psi
  %% h must be only a function of y
  %% tau is constant, along x, usually 1;
  %% along x, we have -x_{max}<=x<=x_{max}
  %% use only first 10 EOF
  xmax=max(x);
  nx=length(x);
  ny=length(h);
  dy=2/(ny-1);
  h_y=firstder(h,dy);
  %% get EOFs Y_n
  ny=length(h);%% expect ny-2 EOF
  nx=length(x);
  dy=2/(ny-1);
				% set up EOF problem: matrices A and B
  A=zeros(ny,ny);
  B=A;
  for n=2:ny-1
    A(n,n-1)=1/((h(n)+h(n-1))/2)^3;
    A(n,n+1)=1/((h(n)+h(n+1))/2)^3;
    A(n,n)=-1/((h(n)+h(n-1))/2)^3-1/((h(n)+h(n+1))/2)^3;
    B(n,n)=1/h(n)^3;
  end
  A=A/dy^2;
  B(1,1)=1/h(1)^3;
  B(ny,ny)=1/h(ny)^3;
				% get eigenvectors and values
  [V,D]=eig(A,B);
  [values,ival]=sort(-diag(D));%% when sorted, the two zeros pop to
  %% the beginning of the list, so chop them out in next two linws
  values=values(3:end);
  Y=V(:,ival(3:end));
		       % normalize so that int Y^2/h^3 has unit amplituede
  for ne=1:ny-2
    factor=0;
    for ky=1:ny
      factor=factor+Y(ky,ne)^2/h(ky)^3*dy;
    end
    Y(:,ne)=Y(:,ne)/sqrt(factor); % normalize
  end
				% compute coefs
 				% compute coefs
  for jx=1:nx
    for ne=1:ny-2
      den=sqrt(values(ne));
      X(jx,ne)=1-cosh(den*x(jx))/cosh(den*max(x));
    end
  end
  for ne=1:ny-2
    a(ne)=0;
    for ky=1:ny
      a(ne)=a(ne)-h_y(ky)*Y(ky,ne)/h(ky)^2/2*dy;
    end
    a(ne)=a(ne)/values(ne);
  end
  psi=zeros(ny,nx);
  for jx=1:nx
    for ne=1:10
      for ky=1:ny
	psi(ky,jx)=psi(ky,jx)+a(ne)*X(jx,ne)*Y(ky,ne);
      end
    end
  end
endfunction

