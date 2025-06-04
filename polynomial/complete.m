function [U,V] = complete(Q,tol) 
% COMPLTE   Complete a nonsquare polynomial matrix to a unimodular matrix
% 
% If Q is a tall polynomial matrix then the command
%   [U,V] = COMPLETE(Q,[TOL])
% produces a unimodular matrix U of the form U = [Q R]. If Q is wide 
% then the unimodular matrix U has the form U = [Q; R]. V is the inverse 
% of U.
% 
% If Q does not have full rank or is not prime then no unimodular
% matrix U exists and an error message follows. Also if Q is square 
% non-unimodular an error is reported.
%
% The optional input argument TOL is the tolerance used for the row or
% column reduction of Q that is part of the algorithm.

%     Author: H. Kwakernaak, August, 1999
%     Copyright 1999 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

% Initialization

if nargin<1,
   error('Not enough input arguments.');
end;
eval('Q = pol(Q);','error(peel(lasterr));');
if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING;
end;
      
[n,k] = size(Q);
wide = 0;
if n == k
   if deg(det(Q)) > 0
      error('Matrix is square non-unimodular.')
   end
elseif k>n
   wide = 1;
   Q = transpose(Q);
   [n,k] = size(Q);
end

% Find the row reduced form

[QQ,r,U,V] = rowred(Q,tol);

if deg(QQ)>0
   error('Matrix cannot be completed to a unimodular matrix.')
elseif rank(QQ(1:k,:))<k
   error('Matrix does not have full rank.')
end

% Redefine U and V

S = [   QQ(1:k,:)   zeros(k,n-k)
      zeros(n-k,k)  eye(n-k,n-k) ];
U = S\U;
V = V*S;

% Construct the completion

C = [Q V(:,k+1:n)];

% Define the output arguments

V = U;
U = C;

if wide
   U = transpose(U);
   V = transpose(V);
end;

%end .. complete
