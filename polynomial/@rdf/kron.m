function C = kron(A,B,tol)
%KRON   Kronecker tensor product of right-den fractions
%
% KRON(A,B) is the Kronecker tensor product of A and B.
% The result is a large matrix formed by taking all possible
% products between the elements of A and those of B.  
%
% KRON(A,B,TOL) works with zeroing specified by the input 
% relative tolerance TOL.
%
% The variable symbols of A,B should be the same. When not, a warning
% is issued and the symbols are changed to the standard one. However, 
% if one symbol is 'z' and the second 'z^-1' then the symbols play
% a role. The resulting symbol is taken form A, no warning being issued.

% See also RDF/TIMES, RDF/MTIMES.

%       Author:  J. Jezek  07-Feb-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $
%                     $ Date 30-May-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

na = nargin;
if na==2,
   tol = PGLOBAL.ZEROING;
elseif na==3,
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
else
   error('Not enough input arguments.');
end;

eval('A = rdf(A); B = rdf(B);', 'error(peel(lasterr));');

[tv,Cv,A,B] = testvf(A,B);
if ~tv, warning('Inconsistent variables.');
end;
[th,Ch,A,B] = testhf(A,B,Cv);
if ~th, warning('Inconsistent sampling periods.');
end;

eval('A = mdf(A,tol);','error(peel(lasterr));');
AB = rdf(kron(A,B.frac.num,tol),tol);
ni = A.s(2);
den = cell(1,ni);
for i = 1:ni,
   den{i} = B.frac.den;
end;
DD = blkdiag(den{1:ni});
C = rdf(AB.frac.num, mtimes(DD,AB.frac.den,tol));

if strcmp(A.p,'prop') & strcmp(B.frac.p,'prop') & ...
      A.frac.tp==B.frac.tp,
   props(C,'prop',A.frac.tp);
end;

if strcmp(PGLOBAL.COPRIME,'cop'), C = coprime(C,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), C = reduce(C,tol);
else C = smreduce(C);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), C = defract(C);
end;

%end .. @rdf/kron
