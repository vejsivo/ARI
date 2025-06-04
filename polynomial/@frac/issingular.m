function test = issingular(A,tol)
%ISSINGULAR  Test if fraction is singular
%
% ISSINGULAR(A) or ISSINGULAR(A,TOL)
% returns 1 if the square fraction A (left-den fraction,
% right-den fraction or scalar=den fraction) is singular
% within tolerance TOL, and 0 otherwise. The tolerance
% TOL is interpreted as for polynomial matrices.
%
% See also: FRAC/RANK, FRAC/ISFULLRANK.

%    Author: J. Jezek  03-Feb-2000
%    Copyright 2000 by Polyx, Ltd.
%    $ Revision 3.0 $  $ Date 04-Apr-2000 $
%                      $ Date 25-Jul-2002 $
%                      $ Date 14-Oct-2002 $

Acl = class(A);
if strcmp(Acl,'frac'),
   error('Invalid 1st argument.');
end;

eval('A = frac(A);','error(peel(lasterr));');
if nargin < 2 | isempty(tol),
   eval('test = issingular(A.num);','error(peel(lasterr));');
else           
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
   eval('test = issingular(A.num,tol);','error(peel(lasterr));');
end;

%end .. @frac/issingular



