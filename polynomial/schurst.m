function [S,U,ind] = schurst(C,tol)
%SCHURST  Ordered Schur decomposition
% 
% The command
%    [S,U,IND] = SCHURST(C) 
% computes the ordered complex Schur decomposition
%    C = USU'
% of the matrix C, such that U is unitary (UU'= I) and S 
% is upper triangular. 
%
% The eigenvalues of C with negative real part appear first along 
% the diagonal of S. The vector IND contains the corresponding column 
% indices of U so that the columns of U(:,IND) span the maximal stable
% invariant subspace of C.
%
% An optional tolerance argument TOL may be specified.
% It is used to determine the sign of the diagonal elements.

%   Author: D. Henrion, December 9, 1998.
%   Modified by D. Henrion, June 2, 1999.
%   Updated to 3.0 by D. Henrion, August 30, 2000.
%   Copyright 1998-2000 by Polyx, Ltd.

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isnumeric(C) | ndims(C)>2,
   error('Invalid 1st argument; must be matrix.');
end;
[rC, cC] = size(C);
if (rC ~= cC),
  error('Matrix must be square.');
end;

[U,T] = schur(C);
[U,T] = rsf2csf(U,T);

if nargin == 1, % default tolerance
  tol = 10.0*eps*max(abs(diag(T)));
end;

% ordering of the Schur decomposition
S = T;
order = 0;
while ~order, % while the order is not correct
 order = 1;
 for k = 1:rC-1,
  if (real(S(k,k)) > real(S(k+1,k+1))),
   order = 0; % diagonal elements to swap
   % complex Givens rotation
   b = S(k,k)-S(k+1,k+1); a = S(k,k+1);
   absa = abs(a);
   if absa == 0, c = 0; s = 1;
   else n = norm([a b]); c = absa/n; s = a/absa*(conj(b)/n);
   end;
   G = [c s;-conj(s) c];
   S(k:k+1,:) = G'*S(k:k+1,:); % row rotation on S
   S(:,k:k+1) = S(:,k:k+1)*G; % column rotation on S
   U(:,k:k+1) = U(:,k:k+1)*G; % column rotation on U
  end; % if
 end; % for k
end; % while

ind = real(diag(S)) < -tol;

%end .. schurst


