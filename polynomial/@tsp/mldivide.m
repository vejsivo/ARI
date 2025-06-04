function R = mldivide(Q,P,var)
%MLDIVIDE(\)   Matrix left divide two-sided polynomials
%
% The command  R = MLDIVIDE(Q,P)  or  R = Q\P  performs left
% matrix division of tsp matrix P by Q. If Q is matrix of numbers
% then the result R is a tsp matrix. If Q is polynomial or tsp
% then the result is a left-denominator fraction in variable 'z'
% or 'z^-1', decided standardly by PGLOBAL.DISCRVAR or by an
% optional input argument VAR, e.g. R = MLDIVIDE(Q,P,'z^-1') .
%
% See also TSP/MRDIVIDE, TSP/LDIVIDE, TSP/RDIVIDE.

%       Author:  J. Jezek  10-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 21-Apr-2000 $
%                     $ Date 22-May-2000 $
%                     $ Date 31-Oct-2000 $
%                     $ Date 25-Jan-2002 $

global PGLOBAL;

if nargin==1,
   error('Not enough input arguments.');
elseif nargin==2 | isempty(var),
   var = PGLOBAL.DISCRVAR;
elseif isa(var,'pol'),
   [vs1,vs2,vd] = size(var);
   if all([vs1,vs2,vd]==1)&(~any(var.c(:,:)-[0,1])),
      var = var.v;
   else
      error('Invalid variable symbol.');
   end;
elseif ~isa(var,'char'),
   error('Invalid variable symbol.');
end;
if ~strcmp(var,'z') & ~strcmp(var,'z^-1') & ~strcmp(var,'zi'),
   error('Invalid variable symbol.');
end;
if strcmp(var,'zi'),
   var = 'z^-1';
end;

eval('P = tsp(P);','error(''Invalid 2nd argument.'');');
if isa(Q,'double') & ndims(Q)==2,
   [s1,s2] = size(Q);
   if s1~=s2,
      error('Denominator matrix is not square.');
   end;
   if ~isfullrank(Q),
      error('Denominator matrix is singular.');
   end;
   Q = inv(Q);
   eval('R = Q*P;','error(peel(lasterr));');
else
   eval ('Q = tsp(Q);','error(''Invalid 1st argument.'');');
   [th,h,P,Q] = testht(P,Q);
   if th==0,
      warning('Inconsistent sampling periods.');
   end;
   if P.o>=Q.o,
      D = Q.p; N = shift(P.p,P.o-Q.o,'z');
   else
      D = shift(Q.p,Q.o-P.o,'z'); N = P.p;
   end;
   D.h = h; N.h = h;
   R = 0;
   eval('R = ldf(D,N);', 'error(peel(lasterr));');
   if strcmp(var,'z^-1'),
      R = reverse(R);
   end;
end;

%end .. @tsp/mldivide

     
      