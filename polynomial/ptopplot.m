function ptopplot(varargin)
%PTOPPLOT  Plots polygonal values sets for a polytope of polynomials
%
% The commmand
%    PTOPPLOT(A0,A1,A2,...,AN,QBOUNDS,OMEGA [,'new'])
% plots the polygonal value set at the generalized (possibly complex)
% frequency OMEGA of the polytope of polynomials
%      A = A0 + Q1*A1 + Q2*A2 + ...+ QN*AN
% A0, A1, A2, ..., AN are given polynomials and Q1, Q2, ..., QN are
% uncertainty intervals [QiMIN,QiMAX] given by
%      QBOUNDS = [ Q1MIN Q1MAX; Q2MIN Q2MAX; ...; QNMIN,QNMAX  ]
% The value set is the set of all possible values that A(Q,OMEGA)
% may assume for the given OMEGA.
%
% If OMEGA is s vector containing a set of complex frequencies then 
% the value sets are evaluated at each frequency.
%
% To open a new figure window include the string 'new' as the last 
% input argument.

%       Author(s):  S. Pejchova 15-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 03-Feb-2000 15:18:34   $
%       $Revision: 3.0 $  $Date: 07-Dec-2000  S. Pejchova   $
%                         $Date: 06-Dec-2002  S. Pejchova   $

ni=nargin; X=[]; omega=0; nw=0;
if ~ni, error('Not enough input arguments.');
elseif ni>1 & ischar(varargin{ni}),
   if strcmp(lower(varargin{ni}),'new'), nw=1; ni=ni-1;
   else, error('Invalid command option.');
   end;
end;   
if ni>3, omega=varargin{ni}; ni=ni-1; end;

omega=omega(:).';
if isempty(omega)|(~isa(omega,'double')),
   error('The set of evaluation frequencies must be a a nonempty vector.');
end;
eval(['X=ptopex(varargin{1:',int2str(ni),'});'],'error(peel(lasterr));');
[rA,cA]=size(varargin{1});
Xomega=polyval(X,omega);
[rX,cA,lo]=size(Xomega);
%
if nw, figure; end;
figure(gcf);
if (~ishold)&(~nw), clf; set(gcf,'name',''); end; 
%
sw = warning; warning off;
for ii=1:rA, for jj=1:cA,
   X_ij=NaN*ones(rX/rA+1,lo);
   for kk=1:lo,
      X_kk=Xomega(ii:rA:rX,jj,kk); cc=[];
      X_kk_s=sort(X_kk); N_kk=norm(X_kk);
      if length(X_kk_s)>1,
         f_X=find(abs(diff(X_kk_s)/N_kk)>eps*1e7);
         f_X=[f_X(:);length(X_kk_s)];
         X_kk_s=X_kk_s(f_X); 
      end;
      eval('cc=convhull(real(X_kk_s),imag(X_kk_s));','c=[];');
      if ~isempty(cc), X_ij(1:length(cc),kk)=X_kk_s(cc); 
      else, X_ij(1:length(X_kk_s),kk)=X_kk_s;
      end;
   end;
   Rmax = ceil(max(max(real(X_ij))));  Rmin = floor(min(min(real(X_ij)))); 
   Imax = ceil(max(max(imag(X_ij))));  Imin = floor(min(min(imag(X_ij)))); 
   if Rmax<=0, Rmax=1; end;
   if Rmin>=0, Rmin=-1; end;
   if Imax<=0, Imax=1; end;
   if Imin>=0, Imin=-1; end;
   
   subplot(rA,cA,(ii-1)*cA+jj);
   plot(real(X_ij),imag(X_ij),[Rmin;Rmax;NaN;0;0],[0;0;NaN;Imin;Imax],':k');
   axis tight;
   if ii==rA & jj==1,
      xlabel('Real Axis');  ylabel('Imag Axis');
   end;
   if rA~=1|cA~=1,
      title(['Entry  (',int2str(ii),',',int2str(jj),')']);
   end;
end; end;
title(['Polytope of polynomials']);
warning(sw); 

%end .. ptopplot
