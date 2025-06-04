function [Q,R] = rdiv(N,D,tol)
%RDIV  Right polynomial matrix division
% The command
%    [Q,R] = RDIV(N,D) 
% computes the polynomial matrix quotient Q and the polynomial 
% matrix remainder R such that
%      N = Q*D + R
% and the degree of R is strictly less than the degree of D.
% Moreover, if D is column reduced then the i-th column degree in R
% is strictly less than the i-th column degree in D. 
% If D is nonsingular and column-reduced then the rational matrix
% R*D^(-1) is strictly proper and the matrices Q and R are unique.
%
% There may be no solution if D is singular. In this case all the 
% entries in Q and R are set to NaN.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: LDIV, COLRED.

%    Author: D. Henrion, November 3, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision 3.0 $  $ Date 29-May-2000  J.Jezek  $
%                      $ Date 18-Jul-2000  J.Jezek  $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin < 2,
 error('Not enough input arguments.');
end;

eval('N = pol(N); D = pol(D);', 'error(peel(lasterr));');

if any(any(isnan(N))) | any(any(isnan(D))),
   error('Polynomial is not finite.');
elseif any(any(isinf(N))) | any(any(isinf(D))),
   error('Polynomial is not finite.');
end;   

[tv,Var,N,D] = testvp(N,D);
if tv==0, warning('Inconsistent variables.');
elseif tv==2, error('Inconsistent variables.');
end;
[th,h,N,D] = testhp(N,D,Var);
if ~th, warning('Inconsistent sampling periods.');
end;

if nargin == 2 | isempty(tol), tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.');
end;

[rN, p] = size(N); [rD, cD] = size(D);

if ~any([rD cD] - [1 1]),      % division by a scalar polynomial
 D = eye(p)*D; cD = p; rD = p; % D as a diagonal matrix
elseif cD ~= p,
   error('Matrices of inconsistent dimensions.');
end;

if (rN == 0) | (p == 0) | (rD == 0),
   %error('Empty polynomial matrices are not allowed.');
   R = N; Q = pol(zeros(rN,rD)); return;  %J.Jezek
end;

if deg(N) < 0,

 % N is the zero matrix
 Q = pol(zeros(rN,rD), Var);
 R = pol(zeros(rN,p), Var);
 solution = 1;
 
elseif deg(D) < 0,

 % D is the zero matrix: division by zero
 warning('Divide by zero.');
 solution = 0;

else

 % when necessary, make D column reduced

 Dlead = lcoef(D, 'col');
 rkDlead = rank(Dlead, tol);

 if rkDlead < p,
  [D,RK,U] = colred(D, tol);
  N = N*U;
 end;

 % build Sylvester resultant matrix

 mui = deg(D, 'col'); mui(mui < 0) = 0; mu = max(mui);
 k = max([0 deg(N, 'col') - mui]);

 A = [eye(p*mu) zeros(p*mu, p*(k+1))];

 indices = zeros(1, p*mu);
 for i = 1:p,
  indices(i+([1:mui(i)]-1)*p) = 1;
 end;
 indices = indices .* (1:p*mu);
 indices = indices(indices > 0);

 A = [A(indices, :); sylv(D, k)];

 Coef = N.Coef;
 B = [Coef(:,:) zeros(rN, p*(k+mu-N.degree))];

 % solve linear system X*A = B using QR decomposition if Sylvester matrix
 % has full row rank and SVD decomposition otherwise
 
 [rA, cA] = size(A);
 R = svd(A);
 rkA = sum(R > max(rA, cA) * max(R) * 1e-8 * tol);

 if rkA == rA, % full row rank

  [Q,R] = qr(A');
  X = B*Q(:,1:rA)/R(1:rA,:)'; % unique solution
  solution = 1;

 else % rank deficient

  [U,R,V] = svd(A);
  R = diag(R); R = diag(1./R(1:rkA));
  X = B*V(:,1:rkA)*R*U(:,1:rkA)';

  % check
  check = norm(X*A-B);
  solution = check < tol;

 end;

 if solution,
  % recover original matrix coefficients

  rX = size(X, 2); sX = length(indices);
  if ~isempty(indices),
   R = zeros(rN, p*mu);
   R(:, indices) = X(:, 1:sX);
   R = pol(R, mu-1);
  else
   R = pol(zeros(rN, p));
  end; 

  Q = pol(X(:, sX+1:rX), k);

  if rkDlead < p, % D was not reduced

   % recover R
   R = xab(U, R, tol);

   solution = ~any(any(isnan(R)));

  end;

  R = pzer(R, tol); pprop(R, Var); R.h = h;
  Q = pzer(Q, tol); pprop(Q, Var); Q.h = h;

 end;

end; % N = 0 ?

if ~solution,
 Q = pol(NaN*ones(rN,rD));
 R = pol(NaN*ones(rN,p));
end;

%end .. rdiv
