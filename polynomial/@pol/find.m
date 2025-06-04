function [I,J,V] = find(X)
%FIND    Find indices of nonzero elements of polynomial vector
%
% The command  I = FIND(X)  returns the indices of the nonzero
% elements in the polynomial vector X.
%
% The command  [I,J] = FIND(X)  returns the row and column indices
% of the nonzero elements in the polynomial matrix X.
%
% The command  [I,J,V] = FIND(X)  also returns a polynomial vector
% containing the nonzero elements in X.
%
% See also FIND (standard MATLAB).

%        Author:  J.Jezek  09-Apr-2001
%        Copyright(c) 2001 by Polyx Ltd.

D = deg(X,'ent');
no = nargout;
if no<=1,
   I = find(D>=0);
else
   [I,J] = find(D>=0);
end;

if no==3,
   L = length(I); DD = deg(X,'mat');
   V = pol(ones(L,DD+1),DD);
   for K = 1:L,
      V.c(K,1,:) = X.c(I(K),J(K),:);
   end;
   V.v = X.v; V.h = X.h;
end;

%end .. @pol/find


