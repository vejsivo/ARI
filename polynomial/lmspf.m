function A = lmispf(B,tol,cofactor)
%LMISPF Polynomial spectral factorization via LMIs
%
% Given a para-Hermitian polynomial matrix B positive definite on
% the stability boundary, the command
%
%    A = LMISPF(B)
%
% computes a stable polynomial matrix A such that B = A'*A.
%
% The spectral factor is computed via LMI optimization.
%
% Note that, in the continuous-time case, B must have no zeros at infinity
% otherwise B is not positive definite on the whole imaginary axis including
% infinity. B has no zeros at infinity for example if B is diagonally reduced
% (see macro DIAGRED).
%
% The relative accuracy for solving LMIs can be specified as an optional
% second input argument. Its default value is the global zeroing tolerance.

% The macro is based on ideas published in H. L. Trentelman, P. Rapisarda
% "New Algorithms for Polynomial J-Spectral Factorization" Mathematics of
% Control, Signals, and Systems, Vol. 12, pp. 24-61, 1999 and also in
% Y. Genin, Y. Nesterov, P. Van Dooren "Positive Transfer Functions and
% Convex Optimization" Proceedings of the European Control Conference,
% Karlsruhe, Germany, 1999.

% Author: Didier Henrion, January 27, 2000.
% Modified by Didier Henrion, February 15, 2000.
% Updated to 3.0 by Didier Henrion, September 1, 2000.
% Modified by Jan Jezek, August 2001, arg checking
% Copyright 2000 Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin < 1,
   error('Not enough input arguments.');
end;

if nargin < 2,
 tol = []; cofactor = [];
elseif nargin < 3, % third argument = spectral co-factorization, see LMISPCOF
 cofactor = [];
end;

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

% Discrete-time TSP format, switch to old format 
if isa(B,'tsp'),
   B = pol(shift(B,deg(B)));
else
   eval('B = pol(B);', 'error(peel(lasterr));');
end;

[rB,cB] = size(B);
if rB~=cB,
   error('Matrix is not square.');
end;

if isempty(B),
   A = pol([]); return;
end
sym = symbol(B); per = B.h;
switch sym,
case {''},
   if isempty(cofactor),
      eval('[A,J] = spf(B,tol);', 'error(peel(lasterr));');
   else
      eval('[A,J] = spcof(B,tol);', 'error(peel(lasterr));');
   end
   if any(diag(J)~=1),
      error('Matrix is not positive definite.');
   end;
   return;
case {'s','p'},
 H = [0 1;1 0];
case {'z^-1','d'},
 H = [1 0;0 -1];
case {'z','q'},
 H = [-1 0;0 1];
end;

if strcmp(sym,'z'), symbol(B,'q');
elseif strcmp(sym,'z^-1'), symbol(B,'d');
end;
if norm(B-B') > tol,
 error('Matrix is not para-Hermitian.');
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

d = deg(B); m = ceil(d/2);
n = size(B,1);

setlmis([]);

varP = lmivar(1,[m*n 1]);

% spectral factorization LMI

if verbose,
 disp('LMISPF: Build LMI.');
end;

LMI = newlmi;

% left hand-side matrix (independent of P)

if H(1,1) == 0, % continuous-time
 Q = zeros((m+1)*n);
 Q(1:n, 1:n) = B{0};
 Q(1:n, 1+m*n:(m+1)*n) = B{m}/2;
 Q(1+m*n:(m+1)*n, 1:n) = B{m}'/2;
 Q(1+m*n:(m+1)*n, 1+m*n:(m+1)*n) = (-1)^m*B{2*m};
 for i = 1:m-1,
  Q(1:n, 1+i*n:(i+1)*n) = B{i}/2;
  Q(1+i*n:(i+1)*n, 1:n) = B{i}'/2;
  Q(1+i*n:(i+1)*n, 1+m*n:(m+1)*n) = (-1)^i*B{m+i}/2;
  Q(1+m*n:(m+1)*n, 1+i*n:(i+1)*n) = (-1)^i*B{m+i}'/2;
 end;
else, % discrete-time
 Q = zeros((m+1)*n);
 Q(1+m*n:(m+1)*n, 1+m*n:(m+1)*n) = B{m};
 for i = 1:m,
  Q(1+m*n:(m+1)*n, 1+(m-i)*n:(m-i+1)*n) = B{m+i}';
  Q(1+(m-i)*n:(m-i+1)*n, 1+m*n:(m+1)*n) = B{m+i};
 end;
end;

% right hand-side matrix (dependent of P)
PP = [zeros(m*n, n) eye(m*n)];
QQ = [eye(m*n) zeros(m*n, n)];
lmiterm([-LMI 1 1 0], Q);
lmiterm([-LMI 1 1 varP], H(1,1)*QQ', QQ);
lmiterm([-LMI 1 1 varP], H(2,1)*QQ', PP, 's');
lmiterm([-LMI 1 1 varP], H(2,2)*PP', PP);

% criterion

LMISYS = getlmis;
nd = decnbr(LMISYS);
c = zeros(1,nd);

for k = 1:nd,
 Pk = defcx(LMISYS,k,varP);
 c(k) = -trace(Pk);
end;

if ~isempty(cofactor),
 c = -c; % spectral co-factorization
end;

% solve LMI

if verbose,
 disp('LMISPF: Solve LMI.');
 disp(['LMISPF: Number of decision variables = ' int2str(nd) '.']);
 disp(['LMISPF: Size of LMI system = ' int2str((m+1)*n) 'x' ...
       int2str((m+1)*n) '.']);
 disp(['LMISPF: Relative accuracy = ' num2str(tol) '.']);
end;

options(1) = tol; % relative accuracy
options(5) = 1; % trace off
[copt,xopt] = mincx(LMISYS,c,options);

if isempty(copt),
 disp('LMISPF: LMI is infeasible !');
 error('Input matrix may not be positive definite on the stability boundary.');
end;

% factor constant (Toeplitz or Hankel) matrix

if verbose,
 disp(['LMISPF: Factor constant matrix of size ' int2str((m+1)*n) ...
  'x' int2str((m+1)*n) '.']);
end;

P = dec2mat(LMISYS,xopt,varP);
L = Q+H(1,1)*QQ'*P*QQ+H(2,1)*PP'*P*QQ+H(1,2)*QQ'*P*PP+H(2,2)*PP'*P*PP;

[U,S,V] = svd(L);
A = pol(sqrt(S(1:n,1:n)) * V(:,1:n)', m, sym);
A.h = per;

if ~isempty(cofactor),
 A = A';
end;

%end .. lmispf

