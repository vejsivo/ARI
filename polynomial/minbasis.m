function [B,rk] = minbasis(A, tol)
%MINBASIS  Minimal polynomial basis
%
% The command
%   [B,RK] = MINBASIS(A) 
% computes a minimal polynomial basis B for the submodule spanned 
% by the columns of the polynomial matrix A.
%
% The basis is minimal in the sense of Forney, that is, B is column 
% reduced and B has full column rank RK. If C is another polynomial 
% basis for A then necessarily deg(C,'col') >= deg(B,'col'). The
% columns of B are arranged according to decreasing degrees.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: COLRED, NULL.

%    Author: D. Henrion, February 9, 1999.
%    Copyright 1998-99 by Polyx, Ltd.

if nargin<1,
   error('Not enough input arguments.');
end;
eval('A = pol(A);', 'error(peel(lasterr));');

if any(any(isnan(A))) | any(any(isinf(A))),
 error('Polynomial is not finite.');
end;

if nargin == 2,
   if ~isa(tol, 'double') | any(size(tol) - [1 1]) | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid 2nd argument.');
   end;
else
   tol = [];
end;

% Extraction of a greatest common divisor of the rows of A

if ~isempty(tol),
 Q = grd(A, tol);
 A = xab(Q, A, tol);
else
 Q = grd(A);
 A = xab(Q, A);
end;

% Reduce to a column-proper matrix

if ~isempty(tol)
 [B,rk] = colred(A,tol);
else
 [B,rk] = colred(A);
end

B = B(:,1:rk);

%end .. minbasis
