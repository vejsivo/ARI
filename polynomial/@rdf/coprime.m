function G = coprime(F,tol);
%COPRIME    Make right-den fraction coprime
%
% Given a right-denominator fraction F, the command
%    G = COPRIME(F)
% cancels a possible right common divisor in the numerator and
% in the denominator of F. The resulting right-den fraction G
% is equal to F but right coprime.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also RDF/REDUCE.

%        Author:  J. Jezek  19-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 26-Apr-2000 $
%                      $ Date 01-Nov-2000 $
%                      $ Date 06-Oct-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

if strcmp(F.frac.c,'cop') & F.frac.tc==tol,    % quick exit
   G = F; return;
end;

Hn = 0; Hd = 0; 
eval('[Hn,Hd] = xayb0(F.frac.den,-F.frac.num,tol);',...
     'error(peel(lasterr));');
[Gn,Gd] = axby0(Hd,-Hn,tol);
G = rdf(Gn,Gd);
G.frac.h = F.frac.h;
props(G,'cop',tol,F.frac.p,F.frac.tp);

%end .. @rdf/coprime
