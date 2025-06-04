function [X,F] = gare(A,B,C,D,E,Q,R,tol);
%GARE    Generalized algebraic Riccati equation
%
% The command 
%    [X,F] = GARE(A,B,C,D,E,Q,R[,TOL]);
% computes the solution X of the generalized algebraic 
% Riccati equation
%    X'*A + A'*X + C'*Q*C - (X'B+C'*D)*R\(B'*X+D'*C) = 0
%    X'*E = E'*X
% and the gain
%    F = R\(B'*X+D'*C)
% such that the feedback law u = -Fx stabilizes the 
% descriptor system 
%    E dx/dt = Ax + Bu
% and makes the impulsive modes non-impulsive. Finite 
% closed-loop poles on the imaginary axis are allowed.
%
% If the descriptor system is not stabilizable or impulse 
% controllable then X is returned as a matrix filled with 
% Infs. In this case F still has a well-defined solution. 
% The corresponding feedback law stabilizes the stabilizable 
% modes and makes the controllable impulsive modes non-impulsive.
%
% The optional tolerance TOL is used by the routine clements 
% and also to test whether the GARE has a finite solution. 
% Its default value is 1e-12.
%
% The algorithm for the solution of the GARE relies on
% transforming the associated Hamiltonian pencil to Clements
% form.

% Huibert Kwakernaak, September-October, 1999
% Copyright 2000 by PolyX Ltd.


% Initialization

if nargin == 7
   tol = 1e-12;
end


% Checks

if nargin < 7
   error('Not enough input arguments.');
end;

if ~isa(A,'double') | ndims(A)>2,
   error('Invalid 1st argument; must be constant matrix.');
end;
[n1,n2] = size(A);
if n1 ~= n2
   error('Inconsistent dimensions of 1st argument.')
else n = n1;
end

if ~isa(B,'double') | ndims(B)>2,
   error('Invalid 2nd argument; must be constant matrix.');
end;
[n1,k] = size(B);
if n1 ~= n
   error('Inconsistent dimensions of 2nd argument.')
end

if ~isa(C,'double') | ndims(C)>2,
   error('Invalid 3rd argument; must be constant matrix.');
end;
[m,n1] = size(C);
if n1 ~= n
   error('Inconsistent dimensions of 3rd argument.')
end

if ~isa(D,'double') | ndims(D)>2,
   error('Invalid 4th argument; must be constant matrix.');
end;
[m1,k1] = size(D);
if m1 ~= m | k1 ~= k
   error('Inconsistent dimensions of 4th argument.')
end

if ~isa(E,'double') | ndims(E)>2,
   error('Invalid 5th argument; must be constant matrix.');
end;
[n1,n2] = size(E);
if n1 ~= n2 | n1 ~= n
   error('Inconsistent dimensions of 5th argument.');
end

if ~isa(Q,'double') | ndims(Q)>2,
   error('Invalid 6th argument; must be constant matrix.');
end;
[m1,m2] = size(Q);
if m1 ~= m2 | m1 ~= m
   error('Inconsistent dimensions of 6th argument.');
end

if ~isa(R,'double') | ndims(R)>2,
   error('Invalid 7th argument; must be constant matrix.');
end;
[k1,k2] = size(R);
if k1 ~= k2 | k1 ~= k
   error('Inconsistent dimensions of 7th argument.');
end

if ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | ...
      tol<0 | tol>1,
   error('Invalid tolerance.');
end;

% Form the Hamiltonian pencil

P12 = s*E-A+B/R*D'*C;
P = [  B/R*B'           P12
        P12'     -C'*Q*C+C'*D/R*D'*C  ];

% Transform to Clements form

[Cl,W,p] = clements1((P+P')/2,0,tol);

% Output

W11 = W(1:n,1:n);
W12 = W(1:n,n+1:2*n);

F = R\(B'*W11'+D'*C*W12');
[U,S,V] = svd(W12');
r = length(find(diag(S)>tol));
Si = zeros(n,n); 
for i = 1:r, Si(i,i) = 1/S(i,i); end
F = F*V*Si*U';

if r < n
   X = inf*ones(n,n);
else
   X = W11'/W12';
end

%end .. gare
