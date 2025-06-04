function test = issingular(A, tol)
%ISSINGULAR  Test if polynomial is singular
%
% ISSINGULAR(A) or ISSINGULAR(A, TOL) returns 1 if the square
% polynomial matrix A is singular within tolerance TOL, and 0 
% otherwise. The tolerance TOL is interpreted as in the macro RANK.
%
% See also: POL/RANK, POL/ISFULLRANK.

%     Author: D. Henrion,  May 10, 1998.
%     Copyright 1998 by Polyx, Ltd.
%     $Revision 3.0$  $Date 11-8-1999  J. Jezek $

eval('A = pol(A);','error(peel(lasterr));');

As = size(A);
if As(1) ~= As(2),
 error('Matrix is not square.');
end;

if nargin < 2 | isempty(tol), % default tolerance
 r = rank(A);
else           % user-supplied tolerance
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
    error('Invalid tolerance.');
 end;
 r = rank(A, tol);
end;

test = (r < As(1));

%end .. @pol/issingular

