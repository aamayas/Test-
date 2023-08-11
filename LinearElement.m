function [C,a,b,c,theta] = LinearElement (x);
  kmax=3;
  theta=atan2(x(2,2)-x(1,2),x(2,1)-x(1,1));
  xieta=[cos(theta),sin(theta);-sin(theta),cos(theta)]*[x(:,1)';x(:,2)'];
  a=xieta(1,2)-xieta(1,3);b=xieta(1,3)-xieta(1,1);c=xieta(2,3)-xieta(2,1);
  canodes=[0,0;1,0;0,1];% nodes of the canonical triangle
  nodes=[1-canodes(:,1)-canodes(:,2),canodes(:,1),canodes(:,2)]*[-b,0;a,0;0,c];
  T=zeros(3,3);T(:,1)=1;T(:,2)=nodes(:,1);T(:,3)=nodes(:,2);
  C=inv(T);
endfunction
