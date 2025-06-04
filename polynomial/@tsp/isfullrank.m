function test = isfullrank(A,tol)
%ISFULLRANK  Test two-sided polynomial is of full rank
%
% ISFULLRANK or ISFULLRANK(A, TOL) returns 1 if the tsp matrix A
% has full row or column rank within tolerance TOL, and 0 otherwise.
% The tolerance TOL is interpreted as in the macro RANK.
%
% See also: TSP/RANK, TSP/ISSINGULAR.

%    Author: J. Jezek  11-8-1999.
%    Copyright 1998 by Polyx, Ltd.

eval('A = tsp(A);','error(peel(lasterr));');
if nargin < 2 | isempty(tol),
   eval('test = isfullrank(A.p);','error(peel(lasterr));');
else           
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
   eval('test = isfullrank(A.p,tol);','error(peel(lasterr));');
end;

%end .. @tsp/isfullrank



