function eta=etaana(tau,x,h);
  %%ID script to compute eta using normal mode expansion
  %% same procedure as in psiana
  %% h must be only a function of y
  %% tau is constant, along x, usually 1;
  %% along x, we have -x_{max}<=x<=x_{max}
  %% use only first 10 EOF
  xmax=max(x);
  nx=length(x);
  ny=length(h);
  dy=2/(ny-1);
  %% get EOFs Y_n
  %% set up A and B
  A=zeros(ny,ny);
  B=A;
  for n=2:ny-1
    A(n,n-1)=((h(n)+h(n-1))/2)^3;
    A(n,n+1)=((h(n)+h(n+1))/2)^3;
    A(n,n)=-(((h(n)+h(n-1))/2)^3+((h(n)+h(n+1))/2)^3);
    B(n,n)=h(n)^3;
  end
  A(1,1:2)=[1,-1];
  A(ny,ny-1:ny)=[-1,1];
  A=A/dy^2;
  B(1,1)=h(1)^3;
  B(ny,ny)=h(ny)^3;
				% get eigenvectors and values
  [V,D]=eig(A,B);
  %% for eta first and last values are rejected
  values=diag(-D);values=values(2:ny-1);Y=V(:,2:ny-1);
  %% normalize so that product has unit amplituede
  for ne=1:ny-2
    factor=0;
    for ky=1:ny
      factor=factor+Y(ky,ne)^2*h(ky)^3*dy;
    end
    Y(:,ne)=Y(:,ne)/sqrt(factor);
  end

  for jx=1:nx
    for ne=1:ny-2
      den=sqrt(values(ne));
      X(jx,ne)=sinh(den*x(jx))/den/cosh(den*max(x));       
    end
  end
  %%[values,X,Y]=etaXY(x,h);
  for ne=1:ny-2
    a(ne)=0;
    for ky=1:ny
      a(ne)=a(ne)+3/2*Y(ky,ne)*h(ky)^2*dy;
    end
  end
  eta=zeros(ny,nx);
  for jx=1:nx
    for ne=1:10
      den=sqrt(values(ne));
      X(jx,ne)=sinh(den*x(jx))/den/cosh(den*max(x));
      for ky=1:ny
	eta(ky,jx)=eta(ky,jx)+a(ne)*X(jx,ne)*Y(ky,ne);
      end
    end
  end
endfunction

