function test = issingular(A,tol)
%ISSINGULAR  Test if matrix-den fraction is singular
%
% ISSINGULAR(A) or ISSINGULAR(A,TOL)
% returns 1 if the square matrix-den fraction A is singular
% within tolerance TOL, and 0 otherwise. The tolerance
% TOL is interpreted as for polynomial matrices.
%
% See also: MDF/RANK, MDF/ISFULLRANK.

%    Author: J. Jezek  03-Feb-2000
%    Copyright 2000 by Polyx, Ltd.
%    $ Revision 3.0 $  $ Date 26-Apr-2000 $
%                      $ Date 06-Nov-2000 $

eval('A = mdf(A);','error(peel(lasterr));');
if nargin < 2 | isempty(tol),
   F = rdf(A);
   eval('test = issingular(F.n);','error(peel(lasterr));');
else           
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
   F = 0;
   eval('F = rdf(A,tol); test = issingular(F.n,tol);', ...
      'error(peel(lasterr));');
end;

%end .. @mdf/issingular
