function [L,U,P] = lu(A,tol)
%LU   LU factorization for polynomial matrices
%
% Given a polynomial matrix A with nonsingular constant term, the command
%    [L,U] = LU(A)
% computes a upper triangular polynomial matrix U whose diagonal entries
% have nonzero real parts, and a unimodular matrix L such that L(0)
% is a "psychologically lower triangular matrix" (i.e., a product of lower
% triangular and permutation matrices), so that
%    A = L*U
% The commmand
%    [L,U,P] = LU(A) 
% returns a unimodular matrix L with lower triangular L(0), an upper 
% triangular matrix U, and a permutation matrix P so that P*A = L*U.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%    Author: D. Henrion, September 25, 1998.
%    Modified by D. Henrion, September 28, 1999.
%    Copyright 1998-99 by Polyx, Ltd.
%    $ Revision $   $ Date 17-Jul-2001  J.Jezek  $
%                   $ arg checking  $

if nargin < 2 | isempty(tol),
   tol = [];
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

if any(any(isnan(A))) | any(any(isinf(A))),
 error('Polynomial is not finite.');
end;

[rA,cA] = size(A); dA = deg(A);

if rA ~= cA,
 error('Matrix is not square.');
elseif rank(polyval(A, 0)) < rA,
 error('Invalid matrix; constant term must be nonsingular.');
end;

% Reduction to triangular form and computation of the inverse
% of the reduction matrix.

[H, V] = tri(A, 'row', tol);
[Vi, detV] = adj(V, tol); detV = pzer(detV);
if deg(detV) > 0,
 error('Reduction to triangular form failed.');
end;
Vi = Vi /polyval(detV,0);

% LU decomposition of constant term.
[L0,U0,P] = lu(polyval(A,0));

L = P*Vi*polyval(V,0)*P'*L0;
U = (U0/polyval(H,0))*H;

%end .. @pol/lu
