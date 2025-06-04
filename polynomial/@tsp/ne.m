function c = ne(A,B,tol)
%NE  (~=) Test if two-sided polynomials unequal
%           A ~= B
%
% C = A~=B does elementwise comparison of the tsp matrices 
% A and B with global zeroing tolerance. A and B must have the same 
% dimensions unless one is a scalar. The scalar is compared 
% with every enrty of the other matrix.
%
% C = NE(A,B) works alike.
% C = NE(A,B,TOL) works with tolerance specified by the input tolerance TOL.
% 
% See also: TSP/EQ.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $ Revision 3.0 $  $ Date 05-Apr-2000 $
%                         $ Date 22-May-2000 $
%                         $ Date 29-Jan-2002 $

global PGLOBAL;

na = nargin;
if na<2, error('Not enough input arguments.'); end;
if na==2, tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

eval('A = tsp(A);','error(''Invalid 1st argument.'');');
eval('B = tsp(B);', ...
   'eval(''B = mdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(B,'mdf'),
   c = 0;
   eval('c = ne(A,B,tol);','error(peel(lasterr));');
   return;
end;   

[th,Ch,A,B] = testht(A,B);
if th==0,
   warning('Inconsistent sampling periods.');
end;

oo = min(A.o,B.o);
AA = shift(A.p,A.o-oo,'z'); BB = shift(B.p,B.o-oo,'z');
eval('c = ne(AA,BB,tol);','error(peel(lasterr));');

%end .. @tsp/ne

   