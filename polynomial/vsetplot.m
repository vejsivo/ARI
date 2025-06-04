function vsetplot(V,varargin)
%VSETPLOT  Plot value set for a parametric polynomial.
%
% The command 
%     VSETPLOT(V [,PLOT_TYPE [,'new']])
% plots the value matrix V computed by macro VSET.
%
% The next input argument PLOT_TYPE can be used to choose the plot type.
%   PLOT_TYPE='lines' - DEFAULT, plots lines.
%   PLOT_TYPE='points' - plots only points.
%
% To open a new figure window include the string 'new' among the 
% input arguments.
%
% See also VSET.

%       Author(s):  S. Pejchova 26-11-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 04-Feb-2000 10:12:34   $
%       $Revision: 3.0 $  $Date: 29-Aug-2000  S. Pejchova   $
%                         $Date: 08-Aug-2001  J. Jezek  $
%                         $Date: 04-Dec-2002  S. Pejchova   $

pl_t=0; nw=0; hw=0;
if nargin<1,
   error('Not enough input arguments.');
end;
if isempty(V) | ~isa(V,'double') | ndims(V)~=2,
   error('The value set matrix must be a nonempty 2-dimensional double.');
end;
if nargin>1,
   for ni=1:nargin-1,
      arg = varargin{ni};
      if ~isa(arg,'char'),
         error(['Invalid ',nth(ni+1),' argument; should be a plot type or ''new'' .']);
      end;
      switch lower(arg),
      case 'lines',
      case 'points',
         pl_t=1;
      case 'new',
         nw=1; 
     otherwise,
         error(['Invalid ',nth(ni+1),' argument; should be a plot type or ''new'' .']);
      end;
   end;
end;

Rmax = max([max(real(V)),0]); Rmin = min([min(real(V)),0]);
Imax = max([max(imag(V)),0]); Imin = min([min(imag(V)),0]);
Rd = 0.05*abs(Rmax-Rmin); Id = 0.05*abs(Imax-Imin);
Rmax=Rmax+Rd;  Rmin=Rmin-Rd;
Imax=Imax+Id;  Imin=Imin-Id;
if ~Rmax, Rmax=1; end;
if ~Imax, Imax=1; end;
if ~Rmin, Rmin=-0.1*(abs(Rmax)); end;
if ~Imin, Imin=-0.1*(abs(Imax)); end;

if nw, figure; end;
figure(gcf);
if ishold, hw=1; 
else, clf, set(gcf,'name','');
end;

map1=prism(size(V,2)+2);
switch pl_t,
case 0,
   plot(real(V),imag(V),[Rmin;Rmax;NaN;0;0],[0;0;NaN;Imin;Imax],':k');
case 1
   if ~hw, axes; end;
   for ii=1:size(V,2),
      p1=patch('XData',real(V(:,ii)),'YData',imag(V(:,ii)),'FaceColor','none','EdgeColor','none');
      set(p1,'Marker','.','MarkerEdgeColor',map1(ii+1,:));
      hold on;
   end;
   plot([Rmin;Rmax;NaN;0;0],[0;0;NaN;Imin;Imax],':k');
end;
axis tight;
xlabel('Real Axis'); ylabel('Imag Axis');
t_l='';
if ~isempty(inputname(1)),  t_l=[' ',inputname(1)]; end;
title(['The Value Set Matrix',t_l,' for a Parametric Polynomial.']);
if ~hw, hold off, end;

%end .. vsetplot
