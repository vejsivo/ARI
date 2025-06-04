function S = lmianalysis(A,varargin)
%LMIANALYSIS  Analysis of stability of a polynomial matrix via LMIs
%
% The command
%   LMIANALYSIS(A)
% returns 1 if polynomial matrix A (square, nonsingular) is stable,
% and 0 otherwise.
%
% Stability is checked via an LMI feasibility problem.
% Complex polynomial matrices (hence complex LMIs) are handled.
% Polynomial matrices with zeros at infinity are not handled. It is thus
% recommended to reduce A to column- or row-reduced form (these forms are
% free of zeros at infinity, see macro COLRED) before calling LMIANALYSIS.
%
% With the syntax LMIANALYSIS(A,'primal'), stability is checked via the primal
% optimization algorithm (default method). With the syntax
% LMIANALYSIS(A,'dual'), stability is checked via the dual feasibility
% algorithm.
%
% An optional tolerance argument can be specified.
%
% The complementary stability region {s=x+i*y:[I I*s']*B*[I I*s] >= 0}
% can be specified with the syntax LMIANALYSIS(A,B).

% Author: Didier Henrion, January 19, 2000.
% Modified by Didier Henrion, February 16, 2000.
%          by Jan Jezek, August 22, 2001,  arg checking
% Copyright 2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin < 1
   error('Not enough input arguments.');
end;
eval('A = pol(A);', 'error(peel(lasterr));');
[n m] = size(A);
d = deg(A);
r = n*d;

% Options

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');
B = []; % stability region
tol = []; % tolerance
primal = 1; % primal or dual problem

invalid = 0;
lv = length(varargin);
if lv > 0
   for i = 1:lv,
      arg = varargin{i};
         if ~isempty(arg),
            if isa(arg, 'char'), % valid strings
               if strcmp(arg, 'dual'),
                  primal = 0;
               elseif strcmp(arg, 'primal'),
                  primal = 1;
               else,
                  error('Invalid command option.');
               end;     
            elseif isa(arg, 'double'), % matrix or scalar
               if ~any(size(arg) - [1 1]), % scalar = tolerance
                  tol = arg;
               else, % matrix argument: stability region
                  B = arg;
               end;
            else,
               error(['Invalid ',nth(i+1),' argument.']);
            end;
         end;
   end;
end;

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

% -------------
% Special cases
% -------------

if n ~= m, % non-square matrix
   error('Matrix is not square.')
elseif issingular(A,tol), % singular matrix
   error('Matrix is singular.')
elseif isempty(d) | d == 0, % empty or constant matrix
   S = logical(1);
   return;
end;

% ------------------------------
% Complementary stability region
% ------------------------------

if isempty(B), % unspecified stability region

 sA = symbol(A);

 if strcmp(sA,'s') | strcmp(sA,'p'),

   % % closed right half-plane /\ large disk centered in 0 of radius RD
   % % i.e., a zero of norm > RD is assumed to be infinite and stable
   % RD = 1e3;
   % B = [0 0 1 0;0 RD^2 0 0;1 0 0 0;0 0 0 -1];

   B = [0 1;1 0];

 elseif strcmp(sA,'z^-1') | strcmp(sA,'d'),

   % closed unit disk

   B = [1 0;0 -1];

 elseif strcmp(sA,'z') | strcmp(sA,'q'),

   % % closed outside of unit disk /\ large disk centered in 0 of radius RD
   % % i.e., a zero of norm > RD is assumed to be infinite and stable
   % RD = 1e3;
   % B = [-1 0 0 0; 0 RD^2 1 0;0 0 0 0; 0 0 0 -1];

   B = [-1 0;0 1];

 end;

end;

k = size(B,1)/2;

B00 = B(1:k,1:k);
B01 = B(1:k,k+1:2*k);
B10 = B(k+1:2*k,1:k);
B11 = B(k+1:2*k,k+1:2*k);

PP = [zeros(r,n) eye(r)];
QQ = [eye(r) zeros(r,n)];
AA = A{:};
NN = null(AA);

if size(NN,2) ~= r,
 error('Matrix is singular.');
end;

sizeLMI = 0;
complexLMI = 0;

if primal == 0,

 % ---------------------------------
 % Build dual LMI (real and complex)
 % ---------------------------------

 if verbose,
 disp('LMIANALYSIS: Build dual LMI system.');
 end;

 setlmis([]);

 % decision variables Pgh

 varP = zeros(k,k);

 for g = 1:k,
  varP(g,g) = lmivar(1,[r 1]); % rxr symmetric matrix
  for h = g+1:k,
   varP(h,g) = lmivar(2,[r r]); % rxr non-symmetric matrix
  end;
 end;

 % constraint  P = P[gh] > 0

 LMI = newlmi;
 for g = 1:k,
  for h = g:k,
   lmiterm([-LMI h g varP(h,g)], 1, 1);
  end;
 end;
 sizeLMI = sizeLMI + r*k;

 % constraint N'*Fd(P)*N < 0
 % this LMI may be complex, so split real and complex parts

 LMI = newlmi;

 % (1,1) real block

 for g = 1:k,
  lmiterm([LMI 1 1 varP(g,g)], real(B00(g,g))*QQ', QQ);
  lmiterm([LMI 1 1 varP(g,g)], real(B10(g,g))*PP', QQ, 's');
  lmiterm([LMI 1 1 varP(g,g)], real(B11(g,g))*PP', PP);
  for h = g+1:k,
   lmiterm([LMI 1 1 varP(h,g)], real(B00(g,h))*QQ', QQ, 's');
   lmiterm([LMI 1 1 varP(h,g)], real(B10(g,h))*PP', QQ, 's');
   lmiterm([LMI 1 1 varP(h,g)], real(B01(g,h))*QQ', PP, 's');
   lmiterm([LMI 1 1 varP(h,g)], real(B11(g,h))*PP', PP, 's');
  end;
 end;

 % is it a complex LMI ?

 if (norm(imag(B)) > 0) | (norm(imag(NN)) > 0),

  % yes, so build another (2,2) real block
  complexLMI = 1;

  for g = 1:k,
   lmiterm([LMI 2 2 varP(g,g)], real(B00(g,g))*QQ', QQ);
   lmiterm([LMI 2 2 varP(g,g)], real(B10(g,g))*PP', QQ, 's');
   lmiterm([LMI 2 2 varP(g,g)], real(B11(g,g))*PP', PP);
   for h = g+1:k,
    lmiterm([LMI 2 2 varP(h,g)], real(B00(g,h))*QQ', QQ, 's');
    lmiterm([LMI 2 2 varP(h,g)], real(B10(g,h))*PP', QQ, 's');
    lmiterm([LMI 2 2 varP(h,g)], real(B01(g,h))*QQ', PP, 's');
    lmiterm([LMI 2 2 varP(h,g)], real(B11(g,h))*PP', PP, 's');
   end;
  end;

  % and a (2,1) imaginary block

  for g = 1:k,
   lmiterm([LMI 2 1 varP(g,g)], imag(B10(g,g))*PP', QQ);
   lmiterm([LMI 2 1 varP(g,g)], imag(B01(g,g))*QQ', PP);
   for h = g+1:k,
    lmiterm([LMI 2 1 varP(h,g)], imag(B00(g,h))*QQ', QQ);
    lmiterm([LMI 2 1 -varP(h,g)], imag(B00(h,g))*QQ', QQ);
    lmiterm([LMI 2 1 varP(h,g)], imag(B10(g,h))*PP', QQ);
    lmiterm([LMI 2 1 -varP(h,g)], imag(B10(h,g))*PP', QQ);
    lmiterm([LMI 2 1 varP(h,g)], imag(B01(g,h))*QQ', PP);
    lmiterm([LMI 2 1 -varP(h,g)], imag(B01(h,g))*QQ', PP);
    lmiterm([LMI 2 1 varP(h,g)], imag(B11(g,h))*PP', PP);
    lmiterm([LMI 2 1 -varP(h,g)], imag(B11(h,g))*PP', PP);
   end;
  end;

  % split outer factor

  lmiterm([LMI 0 0 0], [real(NN) -imag(NN);imag(NN) real(NN)]);
  sizeLMI = sizeLMI + 2*r;

 else % real LMI only

  % standard outer factor

  lmiterm([LMI 0 0 0], NN);
  sizeLMI = sizeLMI + r;

 end;

 % ---------------------------------
 % Solve dual LMI (real and complex)
 % ---------------------------------

 LMISYS = getlmis;

 if verbose,
  if complexLMI,
   disp('LMIANALYSIS: Solve complex dual LMI system.');
  else
   disp('LMIANALYSIS: Solve real dual LMI system.');
  end;
  disp(['LMIANALYSIS: Number of scalar decision variables: ' ...
  int2str(decnbr(LMISYS)) '.']);
  disp(['LMIANALYSIS: Size of LMI system: ' int2str(sizeLMI) 'x' ...
        int2str(sizeLMI) '.']);
 end;

 options(4) = 30; % last iterations
 options(5) = 1; % trace off
 [tmin, xfeas] = feasp(LMISYS, options);

 S = 1;

 if tmin < 0,

  % LMI is feasible
  if verbose,
   disp('LMIANALYSIS: LMI is feasible.');
  end;

 elseif tmin < tol,

  % LMI is marginally infeasible
  if verbose,
   disp('LMIANALYSIS: Cannot establish strict feasibility.');
  end;

 else

  % LMI is infeasible

  if verbose,
   disp('LMIANALYSIS: LMI is infeasible.');
  end;

  S = 0;

 end;

 % -----------------------------
 % Check feasibility of dual LMI
 % -----------------------------

 if S == 1,

  P = zeros(r*k);
  for g = 1:k,
   P(1+(g-1)*r:g*r,1+(g-1)*r:g*r) = dec2mat(LMISYS, xfeas, varP(g,g));
   for h = g+1:k,
    P(1+(h-1)*r:h*r,1+(g-1)*r:g*r) = dec2mat(LMISYS, xfeas, varP(h,g));
    P(1+(g-1)*r:g*r,1+(h-1)*r:h*r) = P(1+(h-1)*r:h*r,1+(g-1)*r:g*r)';
   end;
  end;

  if min(eig(P)) < -tol,
   disp('***********************************************************');
   disp('WARNING: Inconsistency: P is not positive definite');
   disp(['Min eig(P) = ' num2str(min(eig(P)))]);
   disp('***********************************************************');

   S = 0;

  end;

  F = zeros(n+r);
  for g = 1:k,
   for h = 1:k,
    R = P(1+(g-1)*r:g*r,1+(h-1)*r:h*r);
    F = F + B00(h,g)*QQ'*R*QQ + B10(h,g)*PP'*R*QQ + ...
    B10(h,g)'*QQ'*R*PP + B11(h,g)*PP'*R*PP;
   end;
  end;
  F = NN'*F*NN;

  if max(eig(F)) > tol,
   disp('***********************************************************');
   disp('Inconsistency: N''*Fd(P)*N is not negative definite');
   disp(['Max eig(N''*Fd(P)*N) = ' num2str(max(eig(F)))]);
   disp('***********************************************************');

   S = 0;

  end;

  if (S == 1) & verbose,
   disp('LMIANALYSIS: Feasibility check is OK.');
  end;

 end;

else

 AA = A{:};

 setlmis([]);

 if (norm(imag(B)) > 0) | (norm(imag(AA)) > 0)
 
  % -------------------------
  % Build complex primal LMIs
  % -------------------------

  if verbose,
   disp('LMIANALYSIS: Build primal LMI system.');
  end;

  % decision variable X = XR+i*XI

  nd = (n+r-1)*(n+r)/2; var = zeros(nd,1); mat = zeros(n+r);
  for g = 1:n+r, for h = g+1:n+r,
   mat(h,g) = lmivar(1,[1 1]); mat(g,h) = -mat(h,g);
  end; end;
  varXI = lmivar(3,mat); % (n+r)x(n+r) skew-symmetric matrix
  varXR = lmivar(1,[n+r 1]); % (n+r)x(n+r) symmetric matrix

  % constraint  X >= 0

  LMI = newlmi;
  lmiterm([-LMI 1 1 varXR], 1, 1);
  lmiterm([-LMI 2 2 varXR], 1, 1);
  lmiterm([-LMI 2 1 varXI], 1, 1);
  sizeLMI = sizeLMI + 2*(n+r);

  % constraint Trace X(1:n,1:n) = 1
  % as 1+tol > Trace > 1

  LMI = newlmi;
  for g = 1:n,
   e = zeros(n+r,1); e(g) = 1;
   lmiterm([-LMI 1 1 varXR], e', e);
  end;
  lmiterm([LMI 1 1 0], 1);
  sizeLMI = sizeLMI + 1;

  LMI = newlmi;
  for g = 1:n,
   e = zeros(n+r,1); e(g) = 1;
   lmiterm([LMI 1 1 varXR], e', e);
  end;
  lmiterm([-LMI 1 1 0], 1+tol);
  sizeLMI = sizeLMI + 1;

  % constraint F(X) >= 0

  LMI = newlmi;

  % (1,1) real block

  for g = 1:k,
   lmiterm([-LMI g g varXR], real(B00(g,g))*QQ, QQ');
   lmiterm([-LMI g g varXR], real(B10(g,g))*PP, QQ', 's');
   lmiterm([-LMI g g varXI], -imag(B10(g,g))*PP, QQ', 's');
   lmiterm([-LMI g g varXR], real(B11(g,g))*PP, PP');
   for h = g+1:k,
     lmiterm([-LMI h g varXR], real(B00(h,g))*QQ, QQ');
    lmiterm([-LMI h g varXI], -imag(B00(h,g))*QQ, QQ');
    lmiterm([-LMI h g varXR], real(B10(h,g))*PP, QQ');
    lmiterm([-LMI h g varXI], -imag(B10(h,g))*PP, QQ');
    lmiterm([-LMI h g varXR], real(B01(h,g))*QQ, PP');
    lmiterm([-LMI h g varXI], -imag(B01(h,g))*QQ, PP');
    lmiterm([-LMI h g varXR], real(B11(h,g))*PP, PP');
    lmiterm([-LMI h g varXI], -imag(B11(h,g))*PP, PP');
   end;
  end;

  % (2,2) real block

  for g = 1:k,
   lmiterm([-LMI g+k g+k varXR], real(B00(g,g))*QQ, QQ');
   lmiterm([-LMI g+k g+k varXR], real(B10(g,g))*PP, QQ', 's');
   lmiterm([-LMI g+k g+k varXI], -imag(B10(g,g))*PP, QQ', 's');
   lmiterm([-LMI g+k g+k varXR], real(B11(g,g))*PP, PP');
   for h = g+1:k,
    lmiterm([-LMI h+k g+k varXR], real(B00(h,g))*QQ, QQ');
    lmiterm([-LMI h+k g+k varXI], -imag(B00(h,g))*QQ, QQ');
    lmiterm([-LMI h+k g+k varXR], real(B10(h,g))*PP, QQ');
    lmiterm([-LMI h+k g+k varXI], -imag(B10(h,g))*PP, QQ');
    lmiterm([-LMI h+k g+k varXR], real(B01(h,g))*QQ, PP');
    lmiterm([-LMI h+k g+k varXI], -imag(B01(h,g))*QQ, PP');
    lmiterm([-LMI h+k g+k varXR], real(B11(h,g))*PP, PP');
    lmiterm([-LMI h+k g+k varXI], -imag(B11(h,g))*PP, PP');
   end;
  end;

  % (2,1) imaginary block

  for g = 1:k,
   for h = 1:k,
    lmiterm([-LMI k+h g varXR], imag(B00(h,g))*QQ, QQ');
    lmiterm([-LMI k+h g varXI], real(B00(h,g))*QQ, QQ');
    lmiterm([-LMI k+h g varXR], imag(B10(h,g))*PP, QQ');
    lmiterm([-LMI k+h g varXI], real(B10(h,g))*PP, QQ');
    lmiterm([-LMI k+h g varXR], imag(B01(h,g))*QQ, PP');
    lmiterm([-LMI k+h g varXI], real(B01(h,g))*QQ, PP');
    lmiterm([-LMI k+h g varXR], imag(B11(h,g))*PP, PP');
    lmiterm([-LMI k+h g varXI], real(B11(h,g))*PP, PP');
   end;
  end;

  sizeLMI = sizeLMI + 2*k*r;

  % criterion Trace(AA'*AA*X)

  AR = real(AA); AI = imag(AA);

  LMISYS = getlmis;
  nd = decnbr(LMISYS);
  c = zeros(1,nd);

  for g = 1:nd,
   XR = defcx(LMISYS,g,varXR);
   XI = defcx(LMISYS,g,varXI);
   c(g) = trace((AR'*AR+AI'*AI)*XR+(AI'*AR-AR'*AI)*XI);
  end;

  % ------------------------
  % Solve complex primal LMI
  % ------------------------

  if verbose,
   disp('LMIANALYSIS: Solve complex primal LMI system.');
   disp(['LMIANALYSIS: Number of scalar decision variables: ' int2str(nd) '.']);
   disp(['LMIANALYSIS: Size of LMI system: ' int2str(sizeLMI) 'x' ...
        int2str(sizeLMI) '.']);
  end;

  options(1) = tol; % relative accuracy
  options(5) = 1; % trace off
  [copt, xopt] = mincx(LMISYS, c, options);

  if ~isempty(copt), % retrieve complex optimal solution
   XR = dec2mat(LMISYS, xopt, varXR);
   XI = dec2mat(LMISYS, xopt, varXI);
   X = XR + sqrt(-1)*XI;
  end;

 else

  % --------------------
  % Build real primal LMIs
  % --------------------

  if verbose,
   disp('LMIANALYSIS: Build primal LMI system.');
  end;

  % decision variable X

  varX = lmivar(1,[n+r 1]); % (n+r)x(n+r) symmetric matrix

  % constraint  X >= 0

  LMI = newlmi;
  lmiterm([-LMI 1 1 varX], 1, 1);
  sizeLMI = sizeLMI +n+r;

  % constraint Trace X(1:n,1:n) = 1
  % as 1+tol > Trace > 1

  LMI = newlmi;
  for g = 1:n,
   e = zeros(n+r,1); e(g) = 1;
   lmiterm([-LMI 1 1 varX], e', e);
  end;
  lmiterm([LMI 1 1 0], 1);
  sizeLMI = sizeLMI + 1;

  LMI = newlmi;
  for g = 1:n,
   e = zeros(n+r,1); e(g) = 1;
   lmiterm([LMI 1 1 varX], e', e);
  end;
  lmiterm([-LMI 1 1 0], 1+tol);
  sizeLMI = sizeLMI + 1;

  % constraint F(X) >= 0

  LMI = newlmi;

  for g = 1:k,
   lmiterm([-LMI g g varX], real(B00(g,g))*QQ, QQ');
   lmiterm([-LMI g g varX], real(B10(g,g))*PP, QQ', 's');
   lmiterm([-LMI g g varX], real(B11(g,g))*PP, PP');
   for h = g+1:k,
    lmiterm([-LMI h g varX], real(B00(h,g))*QQ, QQ');
    lmiterm([-LMI h g varX], real(B10(h,g))*PP, QQ');
    lmiterm([-LMI h g varX], real(B01(h,g))*QQ, PP');
    lmiterm([-LMI h g varX], real(B11(h,g))*PP, PP');
   end;
  end;
  sizeLMI = sizeLMI + k*r;

  % criterion Trace(AA'*AA*X)

  LMISYS = getlmis;
  nd = decnbr(LMISYS);
  c = zeros(1,nd);

  for g = 1:nd,
   X = defcx(LMISYS,g,varX);
   c(g) = trace(AA'*AA*X);
  end;

  % ---------------------
  % Solve real primal LMI
  % ---------------------

  if verbose,
   disp('LMIANALYSIS: Solve real primal LMI system.');
   disp(['LMIANALYSIS: Number of scalar decision variables: ' int2str(nd) '.']);
   disp(['LMIANALYSIS: Size of LMI system: ' int2str(sizeLMI) 'x' ...
        int2str(sizeLMI) '.']);
  end;

  options(1) = tol; % relative accuracy
  options(5) = 1; % trace off
  [copt, xopt] = mincx(LMISYS, c, options);

  if ~isempty(copt), % retrieve optimal solution
   X = dec2mat(LMISYS, xopt, varX);
  end;

 end;

 % -----------------
 % Check feasibility of primal LMI
 % -----------------

 if isempty(copt),
 
  % LMI is infeasible
  error('Inconsistency: LMI is infeasible !'); 

  elseif copt > tol,

  % criterion is greater than tolerance
  if verbose,
   disp(['LMIANALYSIS: stability: criterion=' num2str(copt) ' > tolerance=' ...
    num2str(tol) '.']);
  end;

  S = logical(1);

 else

  % criterion is less than tolerance
  if verbose,
   disp(['LMIANALYSIS: instability: criterion=' num2str(copt) ' <= tolerance=' ...
    num2str(tol) '.']);
  end;

  S = logical(0);

 end


 if verbose,
  disp(['LMIANALYSIS: X size=' int2str(n+r) ' rank=' int2str(rank(X)) ' cond=' ...
       num2str(cond(X)) '.']);
 end;
 if min(eig(X)) < -tol,
  disp('********************************************************');
  disp('WARNING: Inconsistency: X is not positive semi-definite.');
  disp(['Min eig(X) = ' num2str(min(eig(X)))]);
  disp('********************************************************');
 end;

 F = kron(B00, QQ*X*QQ') + kron(B10, PP*X*QQ') + ...
     kron(B01, QQ*X*PP') + kron(B11, PP*X*PP');
 if min(eig(F)) < -tol,
  disp('***********************************************************');
  disp('WARNING: Inconsistency: F(X) is not positive semi-definite');
  disp(['Min eig(F(X)) = ' num2str(min(eig(F)))]);
  disp('***********************************************************');
 end;

end;

%end .. lmianalysis
