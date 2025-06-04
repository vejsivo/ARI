function zpplot(varargin)
%ZPPLOT  Zero-pole map of a polynomial matrix fraction
%
% The commmand
%    ZPPLOT(F)   
% computes the zeros and poles of the polynomial matrix fraction 
% F and plots them in the complex plane. The zeros are rendered 
% as o's  and the poles are  rendered as x's. 
% If F is a matrix-denominator fraction, then all zeros and poles
% corresponding to each particular entry have a same color.
% The commmand
%    ZPPLOT(N,D) 
% computes and plots the zeros and poles of the polynomial matrix 
% fraction  N/D or D\N, where N - numerator and D - denominator are
% polynomial matrices.
% The commmand
%    ZPPLOT(N)   
% computes the roots (zeros) of the polynomial matrix N and plots 
% them in the complex plane as o's.
%
% To open a new figure window include the string 'new' among the
% input arguments.
%
% To hold the current figure window so that subsequent graphing 
% commands add to the existing graph, include the string 'hold' 
% among the input arguments.
%
% For discrete-time matrices also the unit circle is plotted.

%       Author(s):  S. Pejchova 02-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 07-Dec-1998 15:37:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-2000  S. Pejchova   $
%       $Revision:        $Date: 03-Oct-2002  S. Pejchova   $

% Effect on other properties:
% Output is a standard Matlab plot.

ni=nargin;  tw=0; nw=0; I1=[]; nz=1; P=[]; fw=0; hw=0;
narginchk(1,3);
% error(nargchk(1,3,ni));	%REMOVED IN NEW MATLABS
for ii=1:ni,
   argm=varargin{ii};
   if ~nw & ischar(argm)& strcmp(lower(argm),'new'),
      nw=1; 
  elseif ~hw & ischar(argm)& strcmp(lower(argm),'hold'), 
      hw=1;
   else,
      if isa(argm,'frac'), Z=roots(argm.num); tw=tw+2; fw=1;
      else,  eval('argm=pol(argm);','error(peel(lasterr));'); tw=tw+1; 
      end;
      I1=[I1,strmatch(argm.var,{'z^-1';'d';'z';'q'},'exact')];
      
      switch tw,
      case 1, 
         Z=roots(argm); str_name=['Roots map'];
         if ~isempty(inputname(ii)), 
            str_name=[str_name,' of  ',inputname(ii)]; 
         end;
      case 2,
         str_name=['Zero-pole map'];
         if fw,
            fw = class(argm);
            if strcmp(fw,'mdf'), P=argm; 
            else,  P=roots(argm.den); 
            end;
            str_frac = [' ',upper(fw)]; 
            str_frac=[str_frac(1:end-1),'-'];
            str_frac=strrep(str_frac,'FRA-','');
            str_name=[str_name,' of  ',str_frac,'fraction  ',inputname(ii),];
         else,
            eval('argm=pol(argm);','error(peel(lasterr));')
            eval('frac(varargin{nz},argm);','error(peel(lasterr));')
            P=roots(argm);
            str_name=[str_name,' of fraction '];
            if ~isempty(inputname(nz))& ~isempty(inputname(ii)),
               str_name=[str_name,inputname(nz),'-numerator, ',inputname(ii),'-denominator'];
            end; 
         end;
      case 3,
         error('Too many polynomial inputs.');
      end;
      nz=ii;
   end;
end;
if ~tw, error('Not enough polynomial inputs.'); end;
      
if nw, figure;
elseif ~hw, clf; 
else, hold on,
end;
set(gcf,'name',str_name);
figure(gcf);
if fw & strcmp(fw,'mdf'),
   if isempty(P), Z=[], P=[];
   else,
      Z=zeros(1,0); P1=zeros(1,0);
      for ii=1:(size(P,1)), 
         for jj=1:(size(P,2)),
            Zo=roots(P.num(ii,jj)); Po=roots(P.den(ii,jj));
            if isempty(Zo), 
               Z=[Z,NaN*ones(size(Z,1),1)]; 
            else, 
               Z=[[Z;NaN*ones((size(Zo,1)-size(Z,1)),size(Z,2))],...
                     [Zo;NaN*ones((size(Z,1)-size(Zo,1)),1)]];
            end;
            if isempty(Po), 
               P1=[P1,NaN*ones(size(P1,1),1)]; 
            else, 
               P1=[[P1;NaN*ones((size(Po,1)-size(P1,1)),size(P1,2))],...
                     [Po;NaN*ones((size(P1,1)-size(Po,1)),1)]];
            end;
         end; %for jj=1...
      end; %for ii=1...
      P = P1;
      plot(real(Z),imag(Z),'LineStyle','none','Marker','o','MarkerEdgeColor','auto');
      hold on;
      plot(real(P),imag(P),'LineStyle','none','Marker','*','MarkerEdgeColor','auto');
      if ~hw, hold off; end;
   end;
else,
   plot(real(Z),imag(Z),'bo',real(P),imag(P),'r*');
end;
if isempty(Z)&isempty(P),  axis([-1 1 -1 1]); end
if ~isempty(I1),
   hold on;
   xx=[0:pi/50:2*pi];
   plot(sin(xx),cos(xx),'k:');
   if ~hw, hold off; end;  
   axis equal;
end;
a_x = gca;
Rmax = 1.2*max(abs(get(a_x,'XLim')));
Imax = 1.2*max(abs(get(a_x,'YLim')));
set(a_x,'XLim',[-Rmax Rmax],'YLim',[-Imax Imax]);
line([-Rmax Rmax],[0 0],'Color','k');
line([0 0],[-Imax Imax],'Color','k');
xlabel('Real Axis'); ylabel('Imag Axis');
if tw==2, 
   if isempty(Z)&(~isempty(P)),
      [h_l,o_l]=legend('poles','zeros');   
   else,
      [h_l,o_l]=legend('zeros','poles');
   end;
   if fw & strcmp(fw,'mdf'),
      s1=get(findobj(o_l,'Type','text'),'String');
      y_f=findobj(o_l,'Type','line');
      if length(y_f)>1, 
         if strcmp(s1(1,1),'z') & (get(y_f(1),'Ydata')>get(y_f(2),'Ydata')),
            set(y_f(1),'Marker','o'); set(y_f(2),'Marker','*'); 
         else,
            set(y_f(1),'Marker','*'); set(y_f(2),'Marker','o');
         end;
      end;
      set(y_f,'Color',[0 0 0]);
   end;
end;

%end .. zpplot

      
