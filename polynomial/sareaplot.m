function sareaplot(varargin)
%SAREAPLOT  Plots stability area from the (up to 3-D) array produced by SAREA.
%
%The command 
%       SAREAPLOT(Q1, S [,PLOT_TYPE [,'new']])
% or    SAREAPLOT(Q1,Q2, S [,PLOT_TYPE [,'new']])
% or    SAREAPLOT([Q1,Q2,Q3,] S [,PLOT_TYPE [,'new']])
% plots the stability area given by the (1-D, 2-D or 3-D, respecticely) array S
% produced by macro SAREA. 
% The vectors Q1, Q2, Q3 represent the gride used when creating S by SAREA.
% The S must be length(Q1)-by-length(Q2)-by-length(Q3).
%            or length(Q1)-by-length(Q2)
%            or length(Q1)-by-1
%
% The next input argument PLOT_TYPE can be used to choose the plot type.
%   PLOT_TYPE='surf' - DEFAULT, plots surfaces for each Q3's.
%   PLOT_TYPE='points' - plots only stable points.
%   PLOT_TYPE='both' - plots both types of plots (surfaces and points)
%
% To open a new figure window include the string 'new' as the last 
% input argument.
%
% See also SAREA.

%       Author(s):  S. Pejchova 24-11-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 04-Feb-2000 10:34:34   $
%       $Revision: 2.5 $  $Date: 26-Jul-2000 10:52:34   $
%       $Revision: 3.0 $  $Date: 24-Aug-2000  S. Pejchova   $

ni=nargin;  w_m=0;  pl_t=0; vw2=0; nw=0;
if ni==0, error('Not enough input arguments.');
end;
S = varargin{ni};
%if ni>1 & ischar(S),
   while ischar(S),
      switch lower(S),
      case 'new',
         nw=1; ni=ni-1;
      case 'surf',
         ni=ni-1;
      case 'points',
         pl_t=1; ni=ni-1;
      case 'both',
         pl_t=2; ni=ni-1;
      otherwise,
         error('Invalid string for a plot type or a new figure window.');
      end;
      S = varargin{ni};
   end;
%end;
t_l='';
if ~isempty(inputname(ni)),  t_l=[' ',inputname(ni)]; end;

nq = ni-1;
if isempty(S),
   error('Stability array is empty.');
elseif ~isa(S,'double'),
   error('Invalid type of stability array.');  
elseif nq>3, 
   disp('Not yet ready for multidimensional input.'); return;
end;

[X1,Y1,Z1]=size(S);
n_p=prod([X1,Y1,Z1]);
X1=1:X1; Y1=1:Y1; Z1=1:Z1;
xt=X1; yt=Y1; zt=Z1;
if length(X1)==1, xt=[0,1]; end;
if length(Y1)==1, yt=[0,1]; end;
if length(Z1)==1, zt=[0,1]; vw2=1; end;

x_l='1st parameter'; y_l='2nd parameter'; z_l='3rd parameter';
if nq==3,
   if (length(varargin{3})==length(Z1)),
      Z1=varargin{3};  zt=sort(Z1); vw2=0; 
      if ~isempty(inputname(3)),  z_l=inputname(3); end;
   else,  w_m=1;
   end;
end;
if nq>1,
   if (length(varargin{2})==length(Y1)),
      Y1=varargin{2};  yt=sort(Y1);  
      if ~isempty(inputname(2)),  y_l=inputname(2); end;
   else,  w_m=1;
   end;
end;
if nq,
   if (length(varargin{1})==length(X1)),
      X1=varargin{1};  xt=sort(X1);  
      if ~isempty(inputname(1)),  x_l=inputname(1); end;
   else,  w_m=1;
   end;
end;
if w_m, 
   error('Inconsistent dimensions of vector of parameter values and input stability array.'); 
end;

if nw, figure;
else, clf;
end;
set(gcf,'name','');
figure(gcf);

[Yy1,Xx1,Zz1]=meshgrid(Y1,X1,Z1);
if length(X1)==1|length(Y1)==1, pl_t=1; end;
if pl_t,
   %Plot type - 'points'
    S_st=S(:); S_st(S_st<1)=NaN;
   
    Xst=Xx1(:).*S_st;
    Yst=Yy1(:).*S_st;
    Zst=Zz1(:).*S_st;
    plot3(Xst,Yst,Zst,'.m'); grid on;
    axis([min([0,xt]),max([1,xt]),min([0,yt]),max([1,yt]),min([0,zt]),max([1,zt])]);
    h1=gca;
end; 
if pl_t==2, hold on; end;
if pl_t~=1   
   %Plot type - 'surf'
   Wx1=[]; Wy1=[]; Ws1=[];
   for kk=1:length(Z1),
       Ws=S(:,:,kk);
       Ws(Ws<1)=NaN; 
       Wx=Ws.*Xx1(:,:,kk); Wx1=[Wx1,Wx];
       Wy=Ws.*Yy1(:,:,kk); Wy1=[Wy1,Wy];
       Wz=Ws*Z1(kk);
       Ws1=[Ws1,Wz];
    end;
    surf(Wx1,Wy1,Ws1);
    axis([min([0,xt]),max([1,xt]),min([0,yt]),max([1,yt]),min([0,zt]),max([1,zt])]); 
    h1=gca;
 end;
if pl_t==2, hold off; end; 

if length(xt)<10,
   xtl=mat2str(xt);  xtl=strrep(xtl,' ','|');
   xtl=strrep(xtl,'[',''); xtl=strrep(xtl,']','');
   set(h1,'XTick',xt,'XTickLabel',xtl); 
end;
if length(yt)<10,
   ytl=mat2str(yt);  ytl=strrep(ytl,' ','|');
   ytl=strrep(ytl,'[',''); ytl=strrep(ytl,']','');   
   set(h1,'YTick',yt,'YTickLabel',ytl); 
end;
if length(zt)<10, 
   ztl=mat2str(zt);  ztl=strrep(ztl,' ','|');
   ztl=strrep(ztl,'[',''); ztl=strrep(ztl,']','');
   set(h1,'ZTick',zt,'ZTickLabel',ztl); 
end;
xlabel(x_l); ylabel(y_l); zlabel(z_l);
if vw2,
   view(2); 
   if nq==1, ylabel(' '); end;
else,
   hold on;
   Zz=sum(S,3); Zz=Zz|Zz;
   Z_tick=get(h1,'ZTick');
   if ~Z_tick(1),
      z_par=0.1;
      if length(Z_tick)>1, z_par=0.1*Z_tick(2);  end;
   else,
      z_par=2*Z_tick(1);
   end;
   contour3(Xx1(:,:,1),Yy1(:,:,1),z_par*Zz,1,'-k');
   hold off;
   
end;
title(['Stability Area',t_l]);

%end .. sareaplot
