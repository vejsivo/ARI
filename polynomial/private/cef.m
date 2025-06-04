function [A,irow] = cef(A,tol)
%CEF Column echelon form.
%
%    R = CEF(A) computes a column echelon form of A.
%
%    [R,IROW] = CEF(A) also returns the vector IROW of the indices
%    of the independent rows in A, so that A(IROW,:) is a basis for
%    the range of A.
%
%    The reduction is performed by numerically stable column Householder
%    transformations.
%
%    A tolerance TOL may be specified as an additional input argument.
%    Its default value is the global zeroing tolerance.
%
%    See also: RREF, QR, PRIVATE/NULLREF.

%    Author: D. Henrion,  September 20, 1998.
%    Copyright 1998 by Polyx, Ltd.

global PGLOBAL;

[m,n] = size(A);

if nargin < 2, tol = []; end;
me = max(m,n) * norm(A, 1) * 1e-6;
if isempty(tol),
  tol = me * PGLOBAL.ZEROING;   % default tolerance
else
  tol = me * tol;
end;

j = 1;
dep = 0;
irow = zeros(m,1);
Al = A;

for i = 1:m,

 if ~dep,

  % first row in submatrix
  row = Al(1, :); l = length(row); mu = norm(row);

  if mu > tol, % a non zero row..

   irow(i) = 1; % ..is an independent row

   if l > 1,

    r = row(1);
    if abs(r) < tol, % partial column pivoting
     [void, i1] = max(abs(row)); i1 = i1(1);
     col = Al(:, i1); Al(:, i1) = Al(:, 1); Al(:, 1) = col;
     r = row(i1); row(i1) = row(1);
    end;

    % Householder transformation : zeroes the row but its first component 
    beta = r + sign(r) * mu; v = 1; v(2:l) = row(2:l) / beta;
    beta = -2/(v*v'); w = beta * Al * v';
    A(i:m, j:n) = Al + w * v;
    A(i, j+1:n) = zeros(1, n-j); % exact zeroing

   else dep = 1; % all remaining rows are linearly dependent

   end;

   j = j + 1; % next column

  end;

  if j <= n, Al = A(i+1:m, j:n); end; % new submatrix

 end; % if ~dep

end; 

% indices of independent rows in A
irow = irow .* (1:m)'; irow = irow(irow > 0);

%end .. private/cef
