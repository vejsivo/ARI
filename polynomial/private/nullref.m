function Z = nullref(A,irow)
%NULLREF   Null-space in row echelon form.
%
%    Z = NULLREF(A,IROW) computes a basis in row echelon form for
%    the left null-space of A. Matrix A must be in column echelon form
%    and IROW is a column vector of indices of independent rows in A,
%    as returned by CEF.
%
%    The algorithm performs numerically stable Gaussian eliminations
%    with pivoting.
%
%    See also: PRIVATE/CEF, NULL.

%    Author: D. Henrion, May 19, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    Modified by D. Henrion, August 31, 1998. No tolerance.

[m, p] = size(A);
n = m-length(irow);

% Compute null-space.

Z = zeros(n, m);

if n > 0,

 drow = 1:m; drow(irow) = 0; drow = drow(drow > 0); % dependent row indices

 % solve a series of triangular linear systems
 for i = 1:n,
  j = drow(i); k = irow(irow < j); l = 1:length(k);
  Z(i, [k;j]) = [-A(j, l) / A(k, l) 1];
 end;

end;

%end .. private/nullref
