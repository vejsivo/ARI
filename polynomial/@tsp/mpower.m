function [An,C] = mpower(A,n,tol)
%MPOWER      Matrix power of two-sided polynomial
%                   AN = A^N
%
% AN = A^N or AN = MPOWER(A,N) is A to the N-th power if N is 
% scalar integer and A is square. For nonnegative N, the result
% is also two-sided polynomial. For negative N and nonsingular
% A, the result is tsp, if possible, otherwise scalar-denominator
% fraction. The class of the result can be obtained in an
% optional output argument.
%
% The command works with zeroing activated through the global
% variable PGLOBAL.ZEROING. The command AN = MPOWER(A,N,TOL)
% works with zeroing specified by the input tolerance TOL.
%
% See also TSP/POWER, TSP/MTIMES.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 21-Apr-2000  $
%                         $Date: 31-Oct-2000  $


global PGLOBAL;

na = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else,
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

if ~isa(n,'double') | any(size(n)~=1) | ~isreal(n) | floor(n)~=n,
   error('Invalid power; must be integer scalar.');
end;

if n>=0,
   PP = 0;
   eval('PP  = mpower(A.p,n,tol);','error(peel(lasterr));');
   An = tsp(PP); An.h = A.h;
   An = shift(An,A.o*n);
else
   FF = 0;
   eval('FF = inv(A,tol);','error(peel(lasterr));');
   An = mpower(FF,abs(n),tol);
end;
C = class(An);

%end .. @tsp/mpower
