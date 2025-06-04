function khplot(Amin,varargin)
%KHPLOT  Plots Kharitonov rectangles for interval polynomials
%
% The commmand
%    KHPLOT(Amin,Amax,OMEGA) 
% plots the Kharitonov rectangles for the interval polynomial A
% with interval coefficients A_i=[Amin_i,Amax_i] given by its 
% "bounding polynomials" Amin and Amax, evaluated at the
% frequencies given by the vector OMEGA.
%
% To open a new figure window include the string 'new' as the last 
% input argument.

%       Author(s):  M. Sebek,S. Pejchova 02-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 22-Oct-1998 16:27:34   $
%       $Revision: 3.0 $  $Date: 22-Aug-2000  S. Pejchova   $
%                         $Date: 06-Dec-2002  S. Pejchova   $
%                         $Date: 19-Oct-2009  M. Sebek  $ polyval replaced
%                         by value

ni=nargin; nw=0; omega=0; name_str='';
narginchk(1,4);
% error(nargchk(1,4,ni));	%REMOVED IN NEW MATLABS
Amax=Amin; nam1=inputname(1); nam2=nam1;
if ni>1
   if ischar(varargin{end}),
      if strcmp(lower(varargin{end}),'new'), nw=1; ni=ni-1;
      else, error('Invalid command option.');
      end;
   end;
   switch ni,
   case 1,
   case 2, Amax=varargin{1}; nam2=inputname(2);
   case 3, Amax=varargin{1}; omega=varargin{2}; nam2=inputname(2);
   case 4, error('Invalid input argument.');
   end;
end;
omega=omega(:)';
if isempty(omega)|(~isa(omega,'double')),
   error('The set of evaluation frequencies must be a nonempty vector.');
end;
[rA,cA]=size(Amin);
Amin=Amin.'; Amax=Amax.'; Amin=Amin(:); Amax=Amax(:);
eval('[stab,K1,K2,K3,K4,K5,K6,K7,K8]=kharit(Amin,Amax);', ...
   'error(peel(lasterr));');
Komega=zeros(rA*cA,0);
Kp=[K1,K3,K2,K4,K1];
Kn=[K5,K7,K6,K8,K5];
%
for i=1:length(omega),
   if omega(i)>=0,
%      Komega=[Komega,polyval(Kp,j*omega(i))];
        Komega=[Komega,value(Kp,j*omega(i))]; % MS
   else,
%      Komega=[Komega,polyval(Kn,j*omega(i))];
        Komega=[Komega,value(Kn,j*omega(i))]; % MS
  end
end

if nw, figure; end;
h1=gcf;
figure(h1);
if (~ishold)&(~nw), clf; end; 
%
for ii=1:rA, for jj=1:cA,
   indx=(ii-1)*cA+jj;
   K_ij = Komega(indx,:);    
   Rmax = max(abs(real(K_ij)'));
   Imax = max(abs(imag(K_ij)'));
   if ~Rmax, Rmax=1; end;
   if ~Imax, Imax=1; end;
   K_ij=reshape(K_ij,5,(length(K_ij)/5));
   subplot(rA,cA,indx);
   plot(real(K_ij),imag(K_ij),[-Rmax*1.1;Rmax*1.1;NaN;0;0],[0;0;NaN;-Imax*1.1;Imax*1.1],':k');
   axis tight;
   if ii==rA & jj==1,
      xlabel('Real Axis');  ylabel('Imag Axis');
   end;
   nam3 = [' for an Interval Polynomial'];
   if rA~=1|cA~=1,
      title(['Entry  (',int2str(ii),',',int2str(jj),')']);
      nam3 = [' for Interval Polynomials'];
   end;
end; end;
if ~isempty(nam1)&(~isempty(nam2)), name_str=[nam3,'  [',nam1,',',nam2,']']; end;
name_str=['Kharitonov Rectangles ',name_str];
set(h1,'name',name_str);

%end .. khplot
