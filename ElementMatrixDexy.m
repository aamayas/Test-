function [Dex,Dey] = ElementMatrixDexy (C,a,b,c,theta);
  % Cowper: FEM2d.pdf eqn 37
  [kmax,~]=size(C);
  p=[0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0];p=p(1:kmax);
  q=[0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5];q=q(1:kmax);
  P=repmat(p,kmax,1)+repmat(p',1,kmax);
  Pm1=P-1;Pm1=Pm1(:);Pm1(find(Pm1<0))=0;Pm1=reshape(Pm1,[kmax,kmax]);
  Q=repmat(q,kmax,1)+repmat(q',1,kmax);
  Qm1=Q-1;Qm1=Qm1(:);Qm1(find(Qm1<0))=0;Qm1=reshape(Qm1,[kmax,kmax]);
  Ixi=repmat(p,kmax,1).*(c.^(Q+1).*(a.^(Pm1+1)-(-b).^(Pm1+1)).*factorial(Pm1).*factorial(Q)./factorial(Pm1+Q+2));
  Ieta=repmat(q,kmax,1).*(c.^(Qm1+1).*(a.^(P+1)-(-b).^(P+1)).*factorial(P).*factorial(Qm1)./factorial(P+Qm1+2));
  Dxi=transpose(C)*Ixi*C;
  Deta=transpose(C)*Ieta*C;
  ct=cos(theta);st=sin(theta);
  Dex=Dxi*ct-Deta*st;Dey=Dxi*st+Deta*ct;% depends on element orientation
endfunction
