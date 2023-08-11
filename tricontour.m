## -*- texinfo -*-
## @deftypefn {Function File} {} tricontour (@var{tri}, @var{x}, @var{y}, @var{z}, @var{levels})
## @deftypefnx {Function File} {} tricontour (@var{tri}, @var{x}, @var{y}, @var{z}, @var{levels}, @var{linespec})
## Plot level curves for the values of @code{@var{z}} on a triangular mesh in 2D.
##
## The variable @var{tri} is the triangular meshing of the points
## @code{(@var{x}, @var{y})} which is returned from @code{delaunay}. The
## variable @var{levels} is a vector with the values of the countour levels. If
## @var{levels} is a scalar, then it corresponds to the number of
## level curves to be drawn. If exactly one level curve is desired, list
## the level twice in the vector @var{levels}.
##
## If given, @var{linespec} determines the properties to use for the
## lines.
##
## The output argument @var{h} is the graphic handle to the plot.
## @seealso{plot, trimesh, delaunay}
## @end deftypefn

function h = tricontour (tri, x, y, z, levels, string, varargin)
  if (nargin < 6)
    print_usage ();
  endif

  if isscalar(levels);
    dom=[min(z),max(z)];
    dom=mean(dom)+0.99*(dom-mean(dom));
    levels=linspace(dom(1),dom(2),levels);
  endif

  levels=sort(levels);
  lmin=levels(1);
  lmax=levels(length(levels));

  pData=[];                %% no preallocation
%%  pData=zeros(12000,2);  %% preallocation
  pPoints=0;

  for el=1:length(tri)
    values=[z(tri(el,1)),z(tri(el,2)),z(tri(el,3))];
    minval=min(values);
    maxval=max(values);
    locallevel=levels(minval<=levels+eps); # select the levels to be plotted
    if size(locallevel)>0
      locallevel=locallevel(locallevel<=maxval+eps);
    endif
    for level=locallevel
      points=zeros(1,2);
      npoints=1;
      dl=values-level;
      
      if (abs(dl(1))<=10*eps)
        points(npoints,:)=[x(tri(el,1)),y(tri(el,1))];npoints++;endif
      if (abs(dl(2))<=10*eps)
        points(npoints,:)=[x(tri(el,2)),y(tri(el,2))];npoints++;endif
      if (abs(dl(3))<=10*eps)
        points(npoints,:)=[x(tri(el,3)),y(tri(el,3))];npoints++;endif

      if (npoints<=2)
        if ((dl(1)*dl(2)) < 0)    # intersection between 1st and 2nd point
          points(npoints,:)= ( dl(2)*[x(tri(el,1)),y(tri(el,1))]...
                              -dl(1)*[x(tri(el,2)),y(tri(el,2))])/(dl(2)-dl(1));
          npoints++;
        endif
        if ((dl(1)*dl(3)) < 0)    # intersection between 1st and 3rd point
          points(npoints,:)= ( dl(3)*[x(tri(el,1)),y(tri(el,1))]...
                              -dl(1)*[x(tri(el,3)),y(tri(el,3))])/(dl(3)-dl(1));
          npoints++;
        endif
        if ((dl(3)*dl(2)) < 0)    # intersection between 2nd and 3rd point
          points(npoints,:)= ( dl(2)*[x(tri(el,3)),y(tri(el,3))]...
                              -dl(3)*[x(tri(el,2)),y(tri(el,2))])/(dl(2)-dl(3));
          npoints++;
        endif
      endif
      pData=[pData;points; NaN,NaN ];  %% no preallocation
  %%  pData(pPoints+1:pPoints+npoints,1:2)=[points; NaN,NaN ]; %% preallocation
      pPoints += npoints;
    endfor # level
  endfor # el

  pData=pData(1:pPoints-1,:);

  if (nargout>0)
  size(pData)
    h= plot(pData(:,1),pData(:,2),varargin);
  else
    if length(varargin) == 0      
      plot(pData(:,1),pData(:,2),string,'linewidth',2);
    else
      plot(pData(:,1),pData(:,2),string,'linewidth',varargin{1});
    endif
  endif

endfunction

%!demo
