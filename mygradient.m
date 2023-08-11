function [F_x,F_y]=mygradient(F,dx,dy);
  [ny,nx]=size(F);
  for kx=1:nx
    F_y(:,kx)=firstder(F(:,kx),dy);
  end
  for ky=1:ny
    F_x(ky,:)=firstder(F(ky,:),dx);
  end
endfunction
