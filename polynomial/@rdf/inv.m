function G = inv(F,arg2,arg3)
%INV    Inversion of right-den fraction
%          G = INV(F)
%
% For a right-denominator fraction  F = N * D^-1 ,  the command
% returns  the right-den fraction  G = F^-1  = D * N^-1 .
% The numerator N must be square and nonsingular.

%         Author: J. Jezek  23-Nov-1999
%         Copyright(c) by Polyx, Ltd.
%         $ Revision $  $ Date 25-Apr-2000 $
%                       $ Date 01-Nov-2000 $
%                       $ Date 30-Sep-2002 $
%                       $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING;
if nargin>=2,
   if isa(arg2,'double'),     %possible TOL argument accepted
      tol = arg2;
   elseif ~isa(arg2,'char'),  %possible MET argument ignored
      error('Invalid 2nd argument.');
   end;
end;
if nargin==3,
   if isa(arg3,'double'),
      tol = arg3;
   elseif ~isa(arg3,'char'),
      error('Invalid 3rd argument.');
   end;
end;
if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

G = 0;
eval('G = rdf(F.frac.den,F.frac.num);','error(peel(lasterr));');

if strcmp(PGLOBAL.COPRIME,'cop'), G = coprime(G,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), G = reduce(G,tol);
else G = smreduce(G);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), G = defract(G);
end;

%end .. @rdf/inv
