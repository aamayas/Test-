function derivative=firstder(x,dx);
  derivative=zeros(size(x));
  derivative(2:end-1)=(x(3:end)-x(1:end-2))/2/dx;
  derivative(1)=(-3*x(1)+4*x(2)-x(3))/2/dx;
  derivative(end)=(3*x(end)-4*x(end-1)+x(end-2))/2/dx;
endfunction
