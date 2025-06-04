function test = isfullrank(A,tol)
%ISFULLRANK  Test whether a fraction has full rank.           
%
% ISFULLRANK(A) or ISFULLRANK(A,TOL) returns 1 if the fraction A
% (left-den ftraction, right-den fraction or scalar-den fraction)
% has full rank within tolerance TOL, and 0 otherwise. The rank
% is the same as the rank of the numerator. The tolerance TOL is
% interpreted as for polynomial matrices.
%
% See also: FRAC/RANK, FRAC/ISSINGULAR.

%    Author: J. Jezek  03-Feb-2000.
%    Copyright 2000 by Polyx, Ltd.
%    $ Revision $   $ Date 14-Oct-2002 $

eval('A = frac(A);','error(peel(lasterr));');
if nargin < 2 | isempty(tol),
   eval('test = isfullrank(A.num);','error(peel(lasterr));');
else           
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
   eval('test = isfullrank(A.num,tol);','error(peel(lasterr));');
end;

%end .. @frac/isfullrank



