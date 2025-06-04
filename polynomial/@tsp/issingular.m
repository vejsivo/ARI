function test = issingular(A,tol)
%ISSINGULAR  Test if two-sided polynomial is singular
%
% ISSINGULAR or ISSINGULAR(A, TOL) returns 1 if the square tsp 
% matrix A is singular within tolerance TOL, and 0 otherwise.
% The tolerance TOL is interpreted as in the macro RANK.
%
% See also: TSP/RANK, TSP/ISFULLRANK.

%    Author: J. Jezek  11-8-1999.
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision 3.0 $  $ Date 04-Apr-2000 $

eval('A = tsp(A);','error(peel(lasterr));');
if nargin < 2 | isempty(tol),
   eval('test = issingular(A.p);','error(peel(lasterr));');
else           
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
   eval('test = issingular(A.p,tol);','error(peel(lasterr));');
end;

%end .. @tsp/issingular



