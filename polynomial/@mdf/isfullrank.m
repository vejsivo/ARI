function test = isfullrank(A,tol)
%ISFULLRANK  Test if matrix-den fraction has full rank
%
% ISFULLRANK(A) or ISFULLRANK(A,TOL) returns 1 if the matrix-den
% fraction A has full rank within tolerance TOL, and 0 otherwise.
% The tolerance TOL is interpreted as for polynomial matrices.
%
% See also: MDF/RANK, MDF/ISSINGULAR.

%    Author: J. Jezek  03-Feb-2000.
%    Copyright 2000 by Polyx, Ltd.
%    $ Revision $  $ Date 26-Apr-2000 $
%                  $ Date 06-Nov-2000 $

if nargin==1 | isempty(tol),
   if A.frac.s(1)>=A.frac.s(2),
      F = rdf(A);
   else
      F = ldf(A);
   end;
   test = isfullrank(F.n);
elseif isa(tol,'double'),
   F = 0;
   if A.frac.s(1)>=A.frac.s(2),
      eval('F = rdf(A,tol);','error(peel(lasterr));');
   else
      eval('F = ldf(A,tol);','error(peel(lasterr));');
   end;
   test = isfullrank(F.n,tol);
else
   error('Invalid tolerance.');
end;

%end .. @mdf/isfullrank
