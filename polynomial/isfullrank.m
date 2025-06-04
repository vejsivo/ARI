function test = isfullrank(A, tol)
%ISFULLRANK  Test if polynomial is of full rank
%
% ISFULLRANK(A) or ISFULLRANK(A, TOL) returns 1 if the polynomial 
% matrix A has full row or column rank within tolerance TOL, and 0 
% otherwise. The tolerance TOL is interpreted as in the macro RANK.
%
% See also: POL/RANK, POL/ISSINGULAR.

%    Author: D. Henrion,  May 10, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    $Revision 3.0 $   $Date 11-8-1999    J. Jezek  $

if nargin<1,
   error('Not enough input arguments.');
end;
eval('A = pol(A);','error(peel(lasterr));');
 
As = size(A);

if nargin < 2 | isempty(tol), % default tolerance
  r = rank(A);
else           % user-supplied tolerance;
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
     error('Invalid tolerance.');
  end;
  r = rank(A, tol);
end;

test = (r == min(As));

%end .. isfullrank

