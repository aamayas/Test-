function M = ElementMatrixMe (C,a,b,c)
  %$ Cowper: FEM2d.pdf eqn 35
  [kmax,~]=size(C);
  p=[0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0];p=p(1:kmax);
  q=[0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5];q=q(1:kmax);
  P=repmat(p,kmax,1)+repmat(p',1,kmax);
  Q=repmat(q,kmax,1)+repmat(q',1,kmax);
  I=c.^(Q+1).*(a.^(P+1)-(-b).^(P+1)).*factorial(P).*factorial(Q)./factorial(P+Q+2);
  M=transpose(C)*I*C;
endfunction
