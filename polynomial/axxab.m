function X = axxab(A,B,varargin)
%AXXAB  Symmetric polynomial equation solver
%
% The commmand
%    X = AXXAB(A,B) 
% solves the bilateral symmetric matrix polynomial equation
%    A'X + X'A = B
% where B is a para-Hermitian polynomial matrix, that is,
%    B = B' in continuous-time, and
%    B = BH' + BH in discrete-time.
% The degree offset in the matrix B is evaluated upon cancelation of
% the leading and trailing zero matrix coefficients. If there are no
% zero coefficients then the degree offset is degB/2.
%
% The commmand
%    AXXAB(A,B,'syl') 
% solves the equation with the Sylvester matrix method. This is the 
% default method. It may be used with several modifiers:
%
% * In the continuous-time case 
%      AXXAB(A,B,COLDEG) 
%   computes a solution X of column degrees COLDEG. By default, 
%      COLDEG(i) = MAX(DEG(B)-DEG(A),0)
%   for all column indices i. The macro returns a solution with
%   upper-triangular columnwise leading coefficient matrix.
%
% * In the discrete-time case, the macro computes a solution X of 
%   degree DEG(X) = MAX(DEG(A), DEG(B)). The macro returns a solution with
%   upper-triangular absolute coefficient matrix.
%    
% * By default, the Sylvester matrix algorithm does not perform any
%   preliminary rank check and uses Matlab's built-in sparse linear system
%   solver. The command
%      AXXAB(A,B,'chk')
%   can be used to explicitly perform the preliminary rank computation
%   and use Matlab's built-in standard (non-sparse) linear system solver.
%
% The commmand
%    AXXAB(A,B,'red') 
% solves the equation with polynomial reductions, a version of the Euclidean 
% division algorithm for polynomials. If A is stable and A'^-1 B A^-1 is 
% biproper then the macro computes a solution with upper-triangular columnwise
% leading coefficient matrix (continuous-time) and such that XA^-1 is proper.
%
% If there is no solution of degree as specified above then all the entries 
% in X are set equal to NaN.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: AXYAB.

%    Author: D. Henrion, January 19, 1999.
%    Modified by D. Henrion, February 4, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

global PGLOBAL;

if nargin < 2,
 error('The second input argument is missing.');
end;

A = pol(A);
B = pol(B);

if any(any(isnan(A))) | any(any(isnan(B))),
 error('Invalid Not-a-Number entry in the input matrices.');
elseif any(any(isinf(A))) | any(any(isinf(B))),
 error('Invalid infinite entry in the input matrices.');
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

[rA, cA] = size(A); degA = A.deg;
[rB, cB] = size(B); degB = B.deg;

if (rA ~= cA) | (rB ~= cB) | (rA ~= rB),
 error('Inconsistent dimensions of the input matrices.');
end;

n = rA;

% Data type

C = A+B; typeC = C.var;
typeA = A.var; typeB = B.var;
if (~isempty(typeA) & ~strcmp(typeA, typeC)) |  ...
   (~isempty(typeA) & ~strcmp(typeA, typeC)),
 warning('Inconsistent variables.');
 typeC = PGLOBAL.VARIABLE;
else
 if strcmp(typeC, 's') | strcmp(typeC, 'p'),
  type = 'continuous';
 else
  type = 'discrete';
 end;
end;

% Default options.

tol = [];
coldeg = [];
format = 'tri';
method = 'syl';
check = 'no';

% Handle the case when the macro is called from XAAXB.

if nargin > 2,
 if isa(varargin{1}, 'cell'),
  varargin = varargin{1};
 end;
else
 varargin = [];
end;

% Parse optional input arguments.

invalid = 0;
for i = 1:length(varargin),
  arg = varargin{i};
  if ~isempty(arg),
   if isa(arg, 'char'), % valid strings
    if strcmp(arg, 'syl'),
     method = 'syl';
    elseif strcmp(arg, 'red'),
     method = 'red';
    elseif strcmp(arg, 'tri'),
     format = 'tri';
    elseif strcmp(arg, 'chk'),
     check = 'yes';
    else,
     error('Invalid option.');
    end;     
   elseif isa(arg, 'double'), % matrix or scalar
    if ~any(size(arg) - [1 1]), % scalar = tolerance
     if floor(arg) ~= arg, % non integer = tolerance
      tol = arg;
     else
      coldeg = ones(1, n) * arg; % degree
     end;
    elseif strcmp(type, 'discrete'),
     error('The degree index vector is invalid with discrete-time data.');
    elseif any(sort(size(arg)) - [1 n]),
     error(['The length of the degree index vector must be ' int2str(n) '.']);
    else,
     coldeg = arg;
    end;
   else
    error('Invalid input argument.');
   end;
  end;
end;

% tolerance for zeroing:
% relative to TOL and to elements in A and B

maxA = max(abs(nonzeros(A.coef)));
minB = max(1, min(abs(nonzeros(B.coef))));
if maxA > 0, me = min(minB, 1/maxA);
else me = minB; end;

if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else
 tolzero = tol * me;
end;

% relative tolerance for rank computation:

me = 1e-4;

if isempty(tol),
 tolrank = PGLOBAL.ZEROING * me;
else
 tolrank = tol * me;
end;

% check second input matrix

if ~strcmp(type, 'continuous'),

  % skip zero leading coefficients
  first = 0; % index of the first non-zero matrix coefficient
  while (norm(B{first}) < tolzero) & (first < degB),
   first = first + 1;
  end;

  % skip the first 'first' and last 'first' matrix coefficients
  B = pol(B{first:degB-first}, degB-2*first, symbol(B));
  degB = deg(B);

  % odd degree
  if rem(degB, 2),
   error('The input matrix is not para-Hermitian.');
  end;

end;

if norm(B-B') > tolzero,
 warning('Second input matrix is not para-Hermitian.');
end;

% *********************
% Handle special cases.
% *********************

solution = 1;

if isempty(A) | isempty(B),

 solution = 0; % no solution
 
elseif degB < 0,

 % B is zero, so is X
 X = pol(zeros(n), 0);

elseif degA < 0,

 % A is zero but not B
 solution = 0; % no solution

elseif strcmp(type, 'continuous'),

 % ###################
 % CONTINUOUS-TIME
 % ###################

 if ~isreal(A) | ~isreal(B),
   error('Complex continuous-time algorithms not yet implemented.');
 end;

 if n == 1,

  a = A; b = B; % lowercase variable names
  ra = n; ca = cA; rb = rB; cb = cB;
  dega = degA; degb = degB;

  if strcmp(method, 'red'),
 
   % -----------------------------------------------------
   % CONTINUOUS-TIME
   % A, B SCALAR POLYNOMIALS
   % EUCLIDEAN ALGORITHM
   % -----------------------------------------------------

   a = zeros(1, dega+1); b = zeros(1, degb+1);
   a(:) = A.coef; b(:) = B.coef/2;

   % separation of odd and even parts of a and b

   na = floor(dega/2)+1; nb = floor(degb/2)+1; m = max(na, nb);
   ae = zeros(1,m); ao = ae; be = ae;
   a(1, dega+2) = 0; b(1, degb+2) = 0;

   for i = 1:na,
    ae(i) = a(1, (i-1)*2+1);
    ao(i) = a(1, 2*i);
   end;

   for i = 1:nb,
    be(i) = b(1, (i-1)*2+1); % odd coefs supposed to be zero
   end;

   % forward construction of coefficients

   coefs = []; nit = 1;

   while length(find(ao)) > 0, % while polynomial ao(s) is not zero
    a1 = ao(1); a2 = ae(1);
    if (abs(a1) < tolzero) | (abs(a2) < tolzero),
      error('The first input polynomial is not strictly stable.');
    end;
    coefs = [coefs; a2/a1 be(1)/a2];
    be = [be(2:m) - coefs(nit,2)*ae(2:m) 0];
    aen = ao; ao = [ae(2:m) - coefs(nit,1)*ao(2:m) 0]; ae = aen;
    nit = nit + 1;
   end;

   % backward substitutions

   ue = zeros(1,m); uo = ue;
   ue(1:length(be)) = be/ae(1);

   for i = nit-1:-1:1,
    uon = ue - coefs(i,1)*uo;
    ue = [coefs(i,2) uo(1:m-1)]; uo = uon;
   end;

   % construction of solution x(s)

   x = [];
   for i = 1:m,
    x = [x ue(i) -uo(i)];
   end;

   degx = 2*m - 1;

   X = pzer(pol(x, degx), tolzero);

  else % algorithm

   % ---------------------------------------------------------
   % CONTINUOUS-TIME
   % A, B SCALAR POLYNOMIALS
   % SYLVESTER MATRIX METHOD
   % ---------------------------------------------------------

   if isempty(coldeg),
    degx = max(dega, degb);
   else
    if coldeg < degb - dega,
      error('Invalid expected column degrees.');
    else
      degx = coldeg;
    end;
   end;

   RA = sylv(a, degx);
   b = zeros(1, degb+1); b(:) = B.coef/2;
   b = [b zeros(1,dega+degx-degb)];

   be = []; R = [];
   for i = 1:2:dega+degx+1,
    R = [R RA(:, i)];
    be = [be b(i)];
   end;

   rk = rank(R, tolrank);

   if rk ~= size(R, 2),
    x = be * pinv(R); % least-squares solution
    if norm(x*R-be) > tolrank,
     solution = 0; % no solution
    end;
   else
    x = be / R; % exact solution
   end;
   
   if solution,
    for i = 2:2:degx+1,
     x(i) = -x(i);
    end; 
    X = pzer(pol(x, degx), tolzero);
   end;

  end; % algorithm;

 else % scalar or matrix

  % POLYNOMIAL MATRICES

  % permutation matrix for Kronecker product
  n2 = n^2;
  P = zeros(n2);
  for i = 1:n, for j = i:n,
   k = i+(j-1)*n; l = j+(i-1)*n;
   P(k, l) = 1; P(l, k) = 1;
  end; end;

  % reduced permutation indices
  inde = []; indo = [];
  for i = 1:n, for j = 1:n,
   k = i+(j-1)*n; l = j+(i-1)*n;
   if k >= l,
    if k > l,
      indo = [indo l];
    end;
    inde = [inde l];
   end;
  end; end;

  if strcmp(method, 'red'),

   % ---------------------------------------------------------------
   % CONTINUOUS-TIME
   % POLYNOMIAL MATRICES
   % POLYNOMIAL REDUCTION METHOD
   % ---------------------------------------------------------------

   if verbose,
    disp('AXXAB: Compute polynomial matrix adjoint and determinant.');
   end;

   % transformation into a set of polynomial equations
   [T, detA] = adj(A, tol);
   T = pzer(T, tolzero); detA = pzer(detA, tolzero); degdA = detA.deg;
   B2 = pzer(T'*B*T, tolzero); degB2 = B2.deg;

   % resolution of (n+1)n/2 scalar polynomial equations
   % of the form a* x + y* a = 2b
   X2 = pol(zeros(n, n*(degdA+1)), degdA);
   pprop(X2, typeC);

   if verbose,
    disp('AXXAB: Compute diagonal terms.');
   end;

   % diagonal terms
   for i = 1:n,
    xii = axxab(detA, B2(i, i), 'red', tol);
    % X*inv(A) is not proper
    if xii.deg > degdA,
     error('XA^-1 is not proper.');
    end;
    X2(i, i) = xii;
   end;

   if verbose,
    disp('AXXAB: Compute non-diagonal terms.');
   end;

   % n(n-1) other terms
   if (n > 1),
    for i = 1:n-1,
      for j = i+1:n,
       [xij xji] = axyab(detA, 0, B2(i,j), 0, 'red', tol);
       if max(xji.deg, xij.deg) > degdA,
        error('XA^-1 is not proper.');
       end;
       X2(i,j) = xij; X2(j,i) = xji;
      end;
    end;
   end; % if n

   if verbose,
    residue = norm(detA'*X2+X2'*detA-B2);
    if residue > tolzero,
     warning(['AXXAB: First check. High residue  = ' num2str(residue)]);
    else
     disp('AXXAB: First check is OK.');
    end;
   end;

   if verbose,
    disp('AXXAB: Recovering of original solution.');
   end;

   % recovering of X
   X = xab(T, X2, tol);
   if any(any(isnan(X))),
    error('AXXAB: Recovery of X failed.');
   end;

   if verbose,	
     residue = norm(A'*X+X'*A-B);
     if residue > tolzero,
      warning(['Second check. High residue = ' num2str(residue)]);
     else
      disp('AXXAB: Second check is OK.');
     end;
   end;

   X = pzer(X, tolzero); degX = X.deg;

   if strcmp(format, 'tri'),

     PA = A{degA};
     PX = X{degX};

     % construction of Q such that strictly lower part of Q PA + PX is 0 :
     % resolution of a LSE of order (n-1)n/2 

     if verbose,
      disp('AXXAB: Make column-wise leading coefficient upper-triangular.');
     end;

     ind = [];
     for i = 1:n-1,
       ind = [ind 1+i+(i-1)*n:i*n]; % indices of lower triangular entries
     end;
     M = kron(PA',eye(n)) * (eye(n2) - P);
     M = M(ind, ind); % reduced Kronecker product
     y = -PX(:); y = y(ind); % reduced right hand-side -vecPX

     if norm(y) > tolzero, % if lower triangular entries are not zero

      rk = rank(M, tolrank);

      if rk ~= size(M, 1),
        x = pinv(M) * y;
      else
        x = M \ y; % column vector of elements of Q
      end;

      R = zeros(n2, 1); R(ind) = x;
      Q = zeros(n); Q(:) = R; Q = Q - Q'; % construction of skew-symmetric Q
      X2 = pzer(X + shift(Q, degA-degX) * A, tolzero);

     end;

     % check
     residue = norm(A'*X2+X2'*A-B);
     if residue > tolzero,
      warning('AXXAB: Invalid reduction to triangular form.');
      disp('Column-wise leading coefficient is not triangular.');
     else
      X = X2;
     end;

   end; % triangular

  else % algorithm

   % -------------------------------------------
   % CONTINUOUS-TIME
   % POLYNOMIAL MATRICES
   % SYLVESTER MATRIX METHOD
   % -------------------------------------------

   if isempty(coldeg),
    coldeg = max(degB - degA, 0) * ones(1,n);
   end;
   degX = max(coldeg);
  
   % condition on column degrees
   entB = deg(B,'ent');
   colA = deg(A,'col');
   if max(coldeg) < degB - degA,
     error('Invalid expected column degrees.');
   end; 

   if verbose,
    disp('AXXAB: Build linear system.')
   end;

   X = zeros(n, n*(degX+1));

   % suppression of linearly dependent rows
   % even (1+P) : n^2 rows -> n(n+1)/2 independent rows
   % odd  (1-P) : n^2 rows -> n(n-1)/2 independent rows
   Pe = eye(n2) + P; Po = eye(n2) - P;
   Pe = Pe(inde, :); Po = Po(indo, :); 

   % generator of the extended resultant matrix
   KA = zeros(n2, (degA+2*degX+1)*n2);
   for i = 0:degA,
    KA(:, 1+(degA+degX-i)*n2:(degA+degX-i+1)*n2) = ...
      (-1)^i * kron(eye(n), A{i}');
   end;

   if verbose,
    disp('AXXAB: Build reduced resultant matrix.');
   end;

   % construction of the reduced resultant matrix (left hand-side)
   % and the reduced B vector (right hand-side)
   R = []; Bv = [];
   for i = 0:degA+degX,
     if i <= degB,
       Bi = B{i};
       Bi = Bi(:); % vec(Bi)
     else
       Bi = zeros(n2, 1); % if 2*degA > degB
     end;
     if rem(i,2) % odd
       R = [R; Po * KA(:, 1+(degA+degX-i)*n2:(degA+2*degX-i+1)*n2)];
       Bv = [Bv; Bi(indo, :)];
     else % even
      R = [R; Pe * KA(:, 1+(degA+degX-i)*n2:(degA+2*degX-i+1)*n2)];
        Bv = [Bv; Bi(inde, :)];
     end;
   end; % for i  

   clear KA Pe Po P;

   if strcmp(format, 'tri'),
    % additional constraints to force X(degX) triangular
    Ad = zeros(n*(n-1)/2, n2*(degX+1));
    i = 1;
    for col = 1:n, for row = col+1:n,
      Ad(i, row+(col-1)*n+coldeg(col)*n2) = 1;
      i = i + 1;
    end; end;
    R = [R;Ad]; Bv = [Bv;zeros(n*(n-1)/2,1)];
   end;

   % suppression of zero components corresponding to distinct column degrees
   if any(coldeg - degX),
    active_indices = [];
    i = 0;
    for deg_ = 0:degX,
     for col = 1:n,
      if coldeg(col) >= deg_,
       active_indices = [active_indices [n*i + (1:n)]];
      end;
      i = i + 1;
     end;
    end;
    R = R(:, active_indices);
   else
    active_indices = 1:n2*(degX+1);
   end;

   % resolution of augmented LSE
   sol = zeros(n2*(degX+1),1); 

   if verbose,
    disp(['AXXAB: Solve linear system of size ' int2str(size(R,1)) ...
      'x' int2str(size(R,2)) '.']);
   end;

   if strcmp(check, 'no'),

    if verbose,
     disp('AXXAB: Perform sparse Gaussian elimination.');
    end;

    % Fast computation of a solution
    v = sparse(R) \ Bv;

   else

    if verbose,
     disp('AXXAB: Compute rank of compound matrix.');
    end;

    rk = rank(R, tolrank);

    if rk ~= size(R, 1),

     if verbose,
      disp('AXXAB: Compute least-squares solution.');
     end;

     v = pinv(R) * Bv;

     if norm(R*v - Bv) > tolrank,
      solution = 0; % no solution
     end;

    else

     if verbose,
      disp('AXXAB: Perform Gaussian elimination.');
     end;

     v = R \ Bv;

    end;

   end;

   if solution,
    sol(active_indices) = v;
    X(:) = sol;
    X = pzer(pol(X, degX), tolzero);
   end;

  end; % algorithm

 end; % scalar or matrix

 if solution,

  pprop(X, typeC);

  if verbose,
    residue = norm(A'*X+X'*A-B);
    if residue > tolzero,
     warning(['Final check. High residue = ' num2str(residue)]);
    else
     disp(['AXXAB: Final check is OK.']);
    end;
  end;

 else

  X = pol(NaN*ones(n), 0);
  pprop(X, typeC);

  if verbose,
    disp('AXXAB: No polynomial solution was found.');
  end;

 end;

else % continuous or discrete

 % #################
 % DISCRETE-TIME
 % #################

 if n == 1,

  % SCALAR POLYNOMIALS

  ra = n; ca = cA; rb = rB; cb = cB;
  dega = degA; degb = degB;
  a = zeros(1, degA+1); a(:) = A.coef;
  b = zeros(1, degB+1); b(:) = B.coef;

  degbl = degb / 2;

  if strcmp(method, 'red'),

   % ----------------------
   % DISCRETE-TIME
   % SCALAR POLYNOMIALS
   % ----------------------

   br = [b(degbl+1)/2 b(degbl+2:2*degbl+1)];
   bl = conj([b(degbl+1)/2 b(degbl:-1:1)]);

   degx = max(dega,degbl);
   a = [a zeros(1,degx-dega)]; olda = a;

   % forward construction of coefficients
   backward = []; nit = 0;
   while (dega > 0) | (degbl > 0),
    nit = nit + 1;
    if abs(a(1)) < tolzero,
      error(['The first input argument must be strictly stable.']);
    end;
    if degbl > dega,
      coef = bl(degbl+1)/a(1)';
      for i = 1:degbl,	
        bl(i) = bl(i) - coef*a(degbl+2-i)';
      end;
      bl(degbl+1) = 0;
      backward = [backward; 1 coef degbl];
      degbl = degbl - 1;
    else % degbl <= dega
      coef = a(dega+1)/a(1)';
      for i = 1:dega,
        an(i) = a(i) - coef*a(dega+2-i)';
      end;
      a = [an(1:dega) 0];
      backward = [backward; 3 coef dega];
      dega = dega - 1;
    end;
   end %  while

   % backward substitutions
   x = zeros(1,degx+1); xn = x;
   x(1) = bl(1)/a(1)';
   for i = nit:-1:1,
    type = backward(i,1);
    coef = backward(i,2);
    deg_ = backward(i,3);
    if type == 1,
      x(deg_+1) = x(deg_+1) + coef;
    else % type == 3
      for j = 1:deg_+1,
        xn(j) = x(j) - coef*x(deg_+2-j)';
      end;
      x = xn;
    end; % if type
   end; % for i

   % make constant term real
   if abs(imag(x(1))) > tolzero,
    coef = sqrt(-1)*imag(x(1))/real(olda(1));
    for i = 1:degx+1,
     x(i) = x(i) - coef * olda(i);
    end;
   end;

   X = pzer(pol(x, degx), tolzero);

  else % method

   % ----------------------------------
   % DISCRETE-TIME
   % SCALAR POLYNOMIALS
   % SYLVESTER MATRIX METHOD
   % ----------------------------------

   degx = max(dega, degbl);

   % right hand side matrix
   Bs = [b(degbl+1:degb+1).'; zeros(degx-degbl,1)];

   % reduced resultant matrix
   Q1 = fliplr(a).'; Q2 = a.';
   A1 = zeros(2*degx+1,degx+1);
   A2 = zeros(2*degx+1,degx+1);
   for i = 0:degx,
    A1(1+degx-dega+i:degx+i+1,1+i) = Q1;
    A2(1+degx-i:degx+dega+1-i,1+i) = Q2;
   end;
   A1 = A1(1+degx:1+2*degx,:); % reduction
   A2 = A2(1+degx:1+2*degx,:);

   realflag = isreal(A1) & isreal(Bs);
  
   if realflag, % real matrices

    T = A1+A2;
   
   else % imaginary matrices

    T = [real(A1)+real(A2) imag(A1)+imag(A2);
         -imag(A1)+imag(A2) real(A1)-real(A2)];
    Bs = [real(Bs); imag(Bs)];

    % make constant term real
    select = [1:1+degx 3+degx:2+2*degx];
    T = T(:, select);

   end;

   rk = rank([T Bs], tolrank);

   if rk ~= size(T, 1),
    x = pinv(T) * Bs; % least-squares
    if norm(T*x - Bs) > tolrank,
     solution = 0; % no solution
    end;
   else
    x = T \ Bs; % Gaussian elimination
   end;

   if solution,
    if ~realflag, % recover complex solution
     x = x(1:1+degx)+sqrt(-1)*[0;x(2+degx:1+2*degx)];
    end;
    X = pzer(pol(x.', degx), tolzero);
   end;

  end; % algorithm

 else % scalar or matrix

  % POLYNOMIAL MATRICES

  % permutation matrix for Kronecker product
  n2 = n^2;
  P = zeros(n2);
  for i = 1:n, for j = i:n,
   k = i+(j-1)*n; l = j+(i-1)*n;
   P(k, l) = 1; P(l, k) = 1;
  end; end;

  if strcmp(method, 'red'),

   % ----------------------------
   % DISCRETE-TIME
   % POLYNOMIAL MATRICES
   % POLYNOMIAL REDUCTION METHOD
   % ---------------------------- 

   if ~isreal(A) | ~isreal(B),
    error('Polynomial reduction on complex polynomial matrices not yet implemented.');
   end;

   % transformation into a set of polynomial equations
   if verbose,
    disp('AXXAB: Compute polynomial matrix adjoint and determinant.');
   end;

   [T, detA] = adj(A, tol);
   T = pzer(T, tolzero); detA = pzer(detA, tolzero);
   degT = T.deg;
   detAcjg = detA'; degdetA = detA.deg;
   B2 = T'*B*T; degB2 = 2*degT + degB;

   if verbose,
    disp('AXXAB: Compute diagonal terms.');
   end;

   % resolution of (n+1)n/2 scalar polynomial equations
   % of the form a* x + y* a = bl + br
   degX2 = max(degA, degB2/2); % maximal possible degree
   X2 = pol(ones(n, n*(degX2+1)), degX2);
   pprop(X2, typeC);

   % diagonal terms
   for i = 1:n,
    xii = axxab(detA, B2(i, i), 'red', tol);
    X2(i, i) = xii;
   end;

   if verbose,
    disp('AXXAB: Compute non-diagonal terms.');
   end;

   % n(n-1) other terms
   if (n > 1),
    for i = 1:n-1,
      for j = i+1:n,
        [b1 b2 degb2] = polpart(B2(i, j), degB2/2);
        [xij xji] = axyab(detA, b1, b2, degb2, 'red', tol);
        X2(i, j) = xij;
        X2(j, i) = xji;
      end;
    end;
   end; % if n

   if verbose,
    Ap = detA'; d1 = detA.deg;
    Xp = X2'; d2 = X2.deg;
    d3 = B2.deg / 2;
    m = max([d1 d2 d3]);
    p1 = shift(Ap*X2, m - d1); pprop(p1, typeC);
    p2 = shift(Xp*detA, m - d2); pprop(p2, typeC);
    p3 = shift(B2, m - d3); pprop(p3, typeC);
    residue = norm(p1+p2-p3);
    if residue > tolzero,
     warning(['First check. High residue  = ' num2str(residue)]);
    else
     disp('AXXAB: First check is OK.');
    end;
   end;

   if verbose,
    disp('AXXAB: Recovering of original solution.');
   end;

   % recovering of X
   X = xab(T, X2, tol);
   if any(any(isnan(X))),
    error('The recovery of X failed. Try to modify the tolerance.');
   end;

   if verbose,	
     Ap = A'; d1 = A.deg;
     Xp = X'; d2 = X.deg;
     d3 = B.deg / 2;
     m = max([d1 d2 d3]);
     p1 = shift(Ap*X, m - d1); pprop(p1, typeC);
     p2 = shift(Xp*A, m - d2); pprop(p2, typeC);
     p3 = shift(B, m - d3); pprop(p3, typeC);
     residue = norm(p1+p2-p3);
     if residue > tolzero,
      warning(['Second check. High residue = ' num2str(residue)]);
     else
      disp('AXXAB: Second check is OK.');
     end;
   end;

   X = pzer(X, tolzero); degX = X.deg;

   if strcmp(format, 'tri'),

     if verbose,
      disp('AXXAB: Make absolute coefficient upper-triangular.');
     end;

     % absolute coefficient to be transformed to upper triangular form
     PA = A{0};
     PX = X{0};

     % construction of Q such that strictly lower part of Q PA + PX is 0 :
     % resolution of a LSE of order (n-1)n/2 

     ind = [];
     for i = 1:n-1,
      ind = [ind 1+i+(i-1)*n:i*n]; % indices of lower triangular entries
     end;
     M = kron(PA',eye(n)) * (eye(n2) - P);
     M = M(ind, ind); % reduced Kronecker product
     y = -PX(:); y = y(ind); % reduced right handside -vecPX

     if norm(y) > tolzero, % if lower triangular entries are not zero

      rk = rank([M y], tolrank);

      if rk ~= size(M, 1),
       x = pinv(M) * y; % least-square solution
      else
       x = M \ y; % Gaussian elimination
      end;

      R = zeros(n2, 1); R(ind) = x;
      Q = zeros(n); Q(:) = R; Q = Q - Q'; % construction of skew-symmetric Q
      X = pzer(X + Q*A, tolzero);
     
     end; % norm

   end; % triangular

  else % algorithm

   % ----------------------------------
   % DISCRETE-TIME
   % POLYNOMIAL MATRICES
   % SYLVESTER MATRIX METHOD
   % ----------------------------------

   if verbose,
    disp('AXXAB: Build linear system.');
   end;

   % reduced permutation indices
   inde = []; indo = [];
   for i = 1:n, for j = 1:n,
    k = i+(j-1)*n; l = j+(i-1)*n;
    if k >= l,
      if k > l,
        indo = [indo l];
      end;
      inde = [inde l];
    end;
   end; end;

   % suppression of linearly dependent rows
   Pe = P(inde, :);

   degBr = degB/2;
   degX = max(degA, degBr);

   % right hand side matrix

   V = zeros(n2,1);
   V(:) = B{degBr};
   Bs = Pe*V;
   for i = 1:degX,
    if degBr+i <= degB,
     V(:) = B{degBr+i};
    else
     V = zeros(n2,1);
    end;
    Bs = [Bs;V];
   end;

   if verbose,
    disp('AXXAB: Build reduced resultant matrix.');
   end;

   % reduced resultant matrix
   Ap = A'; degAp = Ap.deg;
   R1 = kron(eye(n),Ap).'; R1 = R1.coef; Q1 = R1(:,:).';
   R2 = kron(eye(n),A)*P; R2 = R2.coef; Q2 = R2(:,:).';
   clear R1 R2 P;

   T1 = zeros((2*degX+1)*n2,(degX+1)*n2);
   T2 = zeros((2*degX+1)*n2,(degX+1)*n2);
   for i = 0:degX,
    T1(1+(degX-degA+i)*n2:(degX-degA+degAp+i+1)*n2,1+i*n2:(i+1)*n2) = Q1;
    T2(1+(degX-i)*n2:(degX+degA+1-i)*n2,1+i*n2:(i+1)*n2) = Q2;
   end;
   clear Q1 Q2;

   T1 = [Pe*T1(1+degX*n2:(degX+1)*n2,:); T1(1+(degX+1)*n2:(2*degX+1)*n2,:)];
   T2 = [Pe*T2(1+degX*n2:(degX+1)*n2,:); T2(1+(degX+1)*n2:(2*degX+1)*n2,:)];
   clear Pe;

   if strcmp(format, 'tri'),
    % n(n-1)/2 additional constraints to force X(0) triangular
    Td = zeros(n*(n-1)/2, (degX+1)*n2);
    for i = 1:n*(n-1)/2,
     Td(i, indo(i)) = 1;
    end;
    T1 = [Td; T1]; T2 = [zeros(n*(n-1)/2, (degX+1)*n2); T2];
    Bs = [zeros(n*(n-1)/2, 1); Bs];
   end;

   realflag = isreal(T1) & isreal(T2) & isreal(Bs);

   if realflag, % real matrices

    T = T1+T2; 

   else % imaginary matrices
 
    T = [real(T1)+real(T2) -imag(T1)+imag(T2);
        imag(T1)+imag(T2) real(T1)-real(T2)];
    Bs = [real(Bs); imag(Bs)];

    % make diagonal of constant term real
    E = eye(n); E = ~E(:)';
    select = [ones(1,n2*(degX+1)) E ones(1,n2*degX)];
    select = select .* (1:2*n2*(degX+1));
    select = select(select > 0);
    T = T(:, select);

   end;
   clear T1 T2;

   % LSE resolution
   X = zeros(n, n*(degX+1));

   if verbose,
    disp(['AXXAB: Solve linear system of size ' int2str(size(T,1)) ...
       'x' int2str(size(T,2)) '.']);
   end;

   if strcmp(check, 'no'),

    if verbose,
     disp('AXXAB: Perform sparse Gaussian elimination.');
    end;

    % Fast computation of a solution
    sol = sparse(T) \ Bs;

   else

    if verbose,
     disp('AXXAB: Compute rank of compound matrix.');
    end;

    rk = rank([T Bs], tolrank);

    if rk ~= size(T, 1),

     if verbose,
      disp('AXXAB: Compute least-squares solution.');
     end;

     sol = pinv(T) * Bs; % least-squares
  
     if norm(T*sol - Bs) > tolrank,
      solution = 0;
     end;

    else

     if verbose,
      disp('AXXAB: Perform Gaussian elimination.');
     end;

     sol = T \ Bs; % Gaussian elimination

    end;

   end;

   if solution,
    if ~realflag, % recover complex solution
     Y = zeros(2*n2*(degX+1), 1);
     Y(select) = sol;
     sol = Y(1:n2*(degX+1)) + sqrt(-1)*Y(1+n2*(degX+1):2*n2*(degX+1));
    end;

    X(:) = sol;
    X = pzer(pol(X, degX), tolzero);

   end;

  end; % algorithm

 end; % scalar or matrix

 if solution,

  pprop(X, typeC);

  if verbose,
   Ap = A'; d1 = A.deg;
   Xp = X'; d2 = X.deg;
   d3 = B.deg / 2;
   m = max([d1 d2 d3]);
   p1 = shift(Ap*X, m - d1); pprop(p1, typeC);
   p2 = shift(Xp*A, m - d2); pprop(p2, typeC);
   p3 = shift(B, m - d3); pprop(p3, typeC);
   residue = norm(p1+p2-p3);
   if residue > tolzero,
    warning(['Final check. High residue = ' num2str(residue)]);
   else
    disp(['AXXAB: Final check is OK.']);
   end;
  end;

 else

  X = pol(NaN*ones(n), 0);
  pprop(X, typeC);

  if verbose,
    disp('AXXAB: No polynomial solution was found.');
  end;

 end;

end; % continuous or discrete



