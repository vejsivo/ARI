function [A, J] = spf(B, varargin)
%SPF  Polynomial spectral factorization
%
% If B is a continuous-time para-Hermitian polynomial matrix,
% that is, B = B', then
%    [A,J] = SPF(B)
% solves the polynomial J-spectral factorization problem, i.e.
%    B = A'JA,  
% where A is Hurwitz (stable). The signature matrix J is diagonal
% of the form 
%        diag(J) = [1, 1, .., 1, -1, -1, .., -1].
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% The commmand
%    SPF(B,'ext') 
% performs factorization with factor extraction by interpolation. This is 
% the default method. A tolerance of TOL*norm(C,1) is used for ordering the 
% complex Schur decomposition of the companion matrix C corresponding to B.
%
% The commmand
%    SPF(B,'are') 
% returns a state-space solution by solving an algebraic Riccati equation. 
% The tolerance TOL is used in the macro AXB for recovering the solution of 
% the ARE.
%
% The commmands
%    [A,J] = SPF(B,'ext','nnc')
%    [A,J] = SPF(B,'are','nnc') 
% return a "nearly non-canonical factorization," that is, the factorization
% is rendered in the form
%     B = A'inv(J)A
% with J diagonal constant but not a signature matrix. It may be singular.
%
% The commmands
%    SPF(B,'syl')
%    SPF(B,'red')
% involve a Newton-Raphson iterative scheme based on the corresponding 
% implementation of function AXXAB, provided B is positive definite on the 
% imaginary axis. The signature matrix J then is assumed to be the identity
% matrix. The iterative scheme is stopped when norm(A'*A-B) is less than 
% TOL*norm(B)*100. The tolerance TOL is also used in the macro AXXAB.
%
% If B is a discrete-time para-Hermitian polynomial matrix, that is,
%    B = BH' + BH
% and B is positive definite on the unit circle, then
%    A = SPF(B)
% solves the polynomial spectral factorization problem, i.e.
%        B = A'A,  A is Schur.
% The factorization is performed with a Newton-Raphson iterative scheme
% based on the macro AXXAB. The solution has its absolute coefficient matrix
% A(0) upper triangular with positive diagonal entries.
%
% The commmand
%    SPF(B,'syl') 
% is based on a Sylvester matrix algorithm. This the default method.
% The commmand
%    SPF(B,'red') 
% is based on polynomial reductions.
%
% The iterative scheme is stopped when norm(A'*A-B) is less than
% TOL*norm(B)*100. The tolerance TOL is also used in the macro AXXAB.
%
% If B is a constant symmetric matrix then A is also constant and is
% obtained together with J from the Schur decomposition of B.

%    Authors: D. Henrion, H. Kwakernaak, S. Pejchova, January 19, 1999.
%    Last modified by D. Henrion, January 8, 2003
%    Copyright 1999-2003 by Polyx, Ltd.

global PGLOBAL;

eval('PGLOBAL.FORMAT;','painit;');
if nargin<1,
   error('Not enough input arguments.');
end;

% Handle the case when GRD is called through GLD with
% transposed arguments.

if nargin > 1,
 if isa(varargin{1}, 'cell'),
  varargin = varargin{1};
 end;
 narg = length(varargin)+1;
else
 narg = 1;
end;

eval('B = pol(B);', 'error(peel(lasterr));');

if any(any(isnan(B))) | any(any(isinf(B))),
 error('Polynomial is not finite.');
end;

typeB = B.var;
if strcmp(typeB, 's') | strcmp(typeB, 'p'),
 type = 'continuous';
elseif strcmp(typeB, 'q') | strcmp(typeB, 'z'),
 type = 'forward';
else
 type = 'backward';
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

[n, cB] = size(B);
degB = deg(B);

if (n ~= cB),
 error('Matrix is not square.');
end;

if n==0,
   A = pol([]); J = []; return;
end;
   
if isinf(degB),
 error('Matrix is singular.');
end;


% Default options.

tol = [];
method = [];
canon = 1;

% Parse optional input arguments.

if narg > 1
 for i = 1:length(varargin),
  arg = varargin{i};
  if ~isempty(arg),
   if isa(arg, 'char'), % valid strings
    if strcmp(arg, 'ext'),
     method = 'ext';
    elseif strcmp(arg, 'are'),
     method = 'are';
    elseif strcmp(arg, 'syl'),
     method = 'syl';
    elseif strcmp(arg, 'red'),
     method = 'red';
    elseif strcmp(arg, 'nnc'),
     canon = 0;
    else,
     error('Invalid command option.');
    end;     
   elseif isa(arg, 'double'), % matrix or scalar
    if ~any(size(arg) - [1 1]), % scalar
     tol = arg; % tolerance
    else, % matrix argument: vector of column indices
     error(['Invalid ',nth(i+1),' argument.']);
    end;
   end;
  end;
 end;
end;

if ~isempty(tol),
   if ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

% tolerances:
% - tolzero: reference absolute tolerance
% - tolfact: relative tolerance for iterative Newton process
% - tolneg: relative tolerance for zero elements in Schur decomposition

% reference absolute tolerance
me = norm(B);
if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else
 tolzero = tol * me;
end;

% relative tolerance for iterative Newton process
tolfact = 1e-4;

% tolerance for neglecting zero elements in Schur decomposition
tolneg = 1e-4;

% Check input matrix.

if ~strcmp(type, 'continuous'),

 % skip zero leading coefficients
 first = 0;
 while (abs(B{first}) < tolzero) & (first < degB),
  first = first + 1;
 end;
 B = shift(B, -first);
 degB = degB - first;

 % odd degree
 if rem(degB, 2),
  error('Matrix is not para-Hermitian.');
 end;

end;

if norm(B-B') > tolzero,
 error('Matrix is not para-Hermitian.');
end;

if degB == 0,

  % --------------------------------------------------------------
  % J-factorization of a constant matrix : B = A'*J*A
  % --------------------------------------------------------------

  if verbose,
   disp('SPF: J-factorization of a constant matrix with Schur decomposition.');
  end;

  B = B{0};
  [A,U] = schur(B); A = A';
  D = diag(U); J = diag(sign(D));
  for i = 1:n,
   if abs(D(i)) > tolzero,
    A(i,:) = A(i,:) * sqrt(abs(D(i)));
   end;
  end;
  [J,I] = sort(diag(-J)); 
  J = -diag(J); A = pol(A(I,:), 0);

elseif strcmp(type, 'continuous'),

 % #################
 % CONTINUOUS-TIME
 % #################

 if isempty(method), method = 'ext'; end;

 J = eye(n);

 if issingular(B, tolzero),

  error('Matrix is singular.')

 else

  % --------------------------------------
  % J-factorization of a polynomial matrix
  % --------------------------------------

  if verbose,
   disp('SPF: Diagonal reduction of B.');
  end;

  % B must be diagonally reduced
  [B, void, Udiagred] = diagred(B);
  degB = deg(B);

  % if diagonally reduced B is a constant matrix,
  % perform factorization of a constant matrix
  % otherwise, perform factorization of a polynomial matrix
  
  if degB < 1,

   [A, J] = spf(B);

  elseif strcmp(method, 'ext'),

   % --------------------------------------------------------
   % factor extraction by interpolation
   % --------------------------------------------------------

   % Bm(s) = D'(-s)B(s)D(s) with non-singular highest coef.
   D = deg(B, 'dia');
   ZL = lcoef(B, 'dia');
   m = max(D); Bm = B;
   for i = 1:n,
     d = m - D(i);
     if d > 0, Bm(i,:) = (-1)^d * shift(Bm(i,:), d); end;
   end;
   for i = 1:n,
     d = m - D(i);
     if d > 0, Bm(:,i) = shift(Bm(:,i), d); end;
   end;

   % Make Bm monic
   B2m = Bm{2*m};
   for i = 0:2*m-1,
    Bm{i} = B2m \ Bm{i};
   end;
   Bm{2*m} = eye(n);

   % block companion matrix
   C = compan(Bm); rC = size(C, 1);

   if verbose,
    disp('SPF: Schur decomposition of companion matrix.');
   end;

   % ordered Complex Schur decomposition C = RUR', R unitary.
   [R,U,ind] = schurst(C);
   diagR = diag(R);
   normR = norm(R);
   
   order = rC/2-n*m+sum(D);

   if sum(ind) < order,

     % some zero eigenvalues are missing, so insert them
     % if there are no purely imaginary eigenvalues
     indzero = (abs(real(diagR))/normR < tolneg) ...
	     & (abs(imag(diagR))/normR >=tolneg);
     if ~any(indzero),
      indzero = abs(diagR)/normR < tolneg;
      ind(1:order) = ind(1:order) | indzero(1:order);
     end;
 
   elseif sum(ind) > order,

     % there may be some additional zero eigenvalues, so skip them
     indzero = abs(diagR)/normR >= tolneg;
     ind = ind & indzero;

  end;

   if sum(ind) ~= order,
    disp('Incorrect dimension of the stable invariant subspace.');
    error('The input matrix may have purely imaginary zeros.');
   end;

   if verbose,
    disp('SPF: Interpolation.');
   end;

   % extraction of eigenvalues with negative real parts
   E = U(:, 1:order);

   % active row indices
   indices = 1:n*(m+1);
   for i = 1:n,
    for k = 0:m-D(i)-1,
     indices = indices(indices ~= i+k*n);
    end;
   end;
   Eh = E(indices, :);

   % interpolation
   Th = null(Eh')';

   % if T(s) is not unique...
   if size(Th,1) ~= n,
    disp('SPF: No full rank factor found when performing extraction.');
    disp('The matrix to be factorized may have purely imaginary zeros');
    disp('or the companion matrix may be ill-conditioned.');
    error('The factorization failed !');
   end; 

   if verbose,
    disp('SPF: Recover extraction factor.');
   end;

   % recovering of extraction factor T(s)
   T2 = zeros(n, n*(m+1)); T = T2;
   T2(:, indices) = Th;
   for i = 1:n,
    for k = 0:D(i),
      T(:, i+k*n) = T2(:, i+(m-D(i)+k)*n); % remove spurious zeros
    end;
   end;
  
   T = pol(real(T), m);
   TL = lcoef(T, 'col');

   if canon, % canonical factorization

    if min(svd(TL)) < 1e-5,
     warning('SPF: Factorization close to noncanonical.');
     disp('Consider nearly non-canonical factorization.');
    end;

    if verbose,
     disp('SPF: Compute canonical factorization.');
    end;

    T = inv(TL) * T;

    [K1,J] = spf(ZL); % constant J-factorization

    A = pzer(K1 * T, tolzero);

   else % "nearly" non-canonical factorization

    if verbose,
     disp('SPF: Compute nearly non-canonical factorization.');
    end;

    [Q,J] = schur((TL/ZL)*TL');
    [J,I] = sort(diag(-J)); 
    J = -diag(J); Q = Q(:,I);

    A = Q' * T;

   end;
  
  elseif strcmp(method, 'are'),

   % ----------------------------------------------------------
   % by solution of an Algebraic Riccati Equation
   % ---------------------------------------------------------- 

   if verbose,
    disp('SPF: Prefactorization.');
   end;

   % prefactorization

   degZ = deg(B, 'dia');
   ZL = lcoef(B, 'dia');

   N = pol(zeros(n), 0);
   k = sqrt(norm(B)) / n;

   for i = 1:n, % construction of stable N(s)
    polyn = diag(fliplr(pascal(degZ(i)+1))); % stable polynomial
    N(i, i) = pol(k*polyn'/polyn(degZ(i)+1), degZ(i)); % correct scaling
   end;   

   % Solve N*(s)M(s) + M*(s)N(s) = B(s)
   % with column degrees of M(s) equal to half diagonal degrees of B(s)

   if verbose,
    disp('SPF: Solve bilateral symmetric equation.');
   end;

   M = axxab(N, B, degZ);
   if any(any(isnan(M))),
     error('The prefactorization failed.');
   end;

   if verbose,
    disp('SPF: Build state-space realization.');
   end;
  
   Q = [N; M];
   delta = deg(Q, 'col');
   E = pol(zeros(n), 0);
 
   for i = 1:n,
    E(i, i) = pol([zeros(1,delta(i)) 1], delta(i));
   end;

   if verbose,
    disp('SPF: Compute minimal realization.');
   end;

   % minimal realization of H=Q/E
   [As,Bs,Cs,Ds] = rmf2ss(Q, E);

   if verbose,
    disp('SPF: Schur factorization of Hamiltonian.');
   end;

   % Hamiltonian
   W = [zeros(n) eye(n); eye(n) zeros(n)];
   L = Ds'*W*Ds;
   if min(svd(L)) < tolzero,
     disp('The term L=D''WD in the ARE is singular.');
     error('The input matrix may be singular.');
   end;

   Li = inv(L); S = Cs'*W*Ds;
   H = [As-Bs*Li*S' -Bs*Li*Bs'; -Cs'*W*Cs+S*Li*S' -(As-Bs*Li*S')'];
   [R,U,ind] = schurst(H);
   normR = norm(R);
   
   order = size(As,1);

   if sum(ind) < order,

     % some zero eigenvalues are missing, so insert them
     % if there are no purely imaginary eigenvalues

     indzero = (abs(real(diag(R)))/normR <tolneg) ...
	     & (abs(imag(diag(R)))/normR >=tolneg);
     if ~any(indzero),
      indzero = abs(diag(R))/normR < tolneg;
      ind(1:order) = ind(1:order) | indzero(1:order);
     end;

   elseif sum(ind) > order,

     % there may be some additional zero eigenvalues, so skip them
     indzero = abs(diag(R))/normR >= tolneg;
     ind = ind & indzero;

   end;

   if sum(ind) ~= order,
    disp('Incorrect dimension of stable invariant subspace of Hamiltonian.');
    error('The input matrix may have purely imaginary zeros.');
   end;

   if verbose,
    disp('SPF: Recover solution of ARE.');
   end;

   % solution of ARE

   X1 = U(1:order,ind); X2 = U(order+1:2*order,ind);

   R = pol([zeros(order) eye(order)], 1) - As; % R(s) = sI-A
   K = axb(R, Bs*E, tol); % (sI-A)K(s)=BE(s)

   if canon, % canonical factorization

    if min(svd(X1)) < 1e-5,
     warning('SPF: X1 is close to singular.');
     disp('Consider nearly non-canonical factorization.');
    end;

    if verbose,
     disp('SPF: Canonical factorization.');
    end;

    X = real(X2/X1);
    F = Li*(Bs'*X+S');
    T = E + F*K;

    [L,J] = spf(L); % constant J-factorization

    A = L*T; 
 
   else % "nearly" non-canonical factorization

    if verbose,
     disp('SPF: Nearly non-canonical factorization.');
    end;

    N = real(null([Bs'*X2+S'*X1;X1]')');
    m = size(Bs,2);
    X1h = N(:,1:m); Gh = -N(:,m+1:m+order); 
 
    T = X1h*L*E + Gh*K;
    
    TL = lcoef(T, 'col');

    [Q,J] = schur((TL/ZL)*TL');
    [J,I] = sort(diag(-J)); 
    J = -diag(J); Q = Q(:,I);

    A = Q'*T;

   end;

  elseif strcmp(method, 'syl') | strcmp(method, 'red'),

   % ---------------------------------------------------------------
   % Newton-Raphson iterative scheme with AXXAB
   % ---------------------------------------------------------------

   % initial matrix A0

   if verbose,
    disp('SPF: Check positiveness of input matrix.');
   end;

   Bs0 = B{0};
   % halves of diagonal degrees of B(s) and leading coefficient matrix
   Bsh = lcoef(B, 'dia');
   p = deg(diag(B), 'row') / 2; 

   if (min(eig(Bs0)) < -tolzero) | (min(eig(Bsh)) < -tolzero),

    error('Matrix is not positive-definite on the imaginary axis.');

   else

    e0 = min(eig(Bs0)); eh = min(eig(Bsh));
    if (min([e0 eh]) <= tolzero) & strcmp(method, 'red'),
     error('Matrix is not positive-definite on the imaginary axis.');
    end; 

    if e0 <= tolzero, % Bs0 only semi-definite positive

     % Schur form should be diagonal : positive eigenvalues first
     [S,Psi,ind] = schurst(-Bs0); S = -S;
     Psi(:,ind) = Psi(:,ind) * diag(sqrt(diag(S(ind,ind))));
     Psi(:,sum(ind)+1:n) = zeros(n,n-sum(ind));

    else

     Psi = chol(Bs0);

    end;

    if eh <= tolzero, % Bsh only semi-definite positive

     % Schur form should be diagonal : positive eigenvalues first
     [S,Phi,ind] = schurst(-Bs0); S = -S;
     Phi(:,ind) = Phi(:,ind) * diag(sqrt(diag(S(ind,ind))));
     Phi(:,sum(ind)+1:n) = zeros(n,n-sum(ind));
    
    else
        
     Phi = chol(Bsh);

    end;

   end;

   if verbose,
     disp('SPF: Build initial factor.');
   end;

   degA = max(p);
   newA = pol(zeros(n), 0);
 
   for i = 1:n,
    for j = 1:n,
     if (i == j), % diagonal terms
      % coefficients of the polynomial (a+bs)^p(i)
      % aii(s) = (psi(i,i)^(1/p(i))+phi(i,i)^(1/p(i))*s)^p(i)
      if abs(Phi(i,i)) < tolzero,
       newA(i,i) = Psi(i,i);
      else
       root = (Psi(i,i)/Phi(i,i))^(1/p(i));
       coef = fliplr(poly(-ones(1,p(i))*root));
       newA(i,i) = pol(coef, p(i));
      end;
     else % non-diagonal terms
      a = [zeros(1,p(j)) Phi(i,j)];
      a(1) = a(1) + Psi(i,j);
      % aij(s) = psi(i,j) + phi(i,j)*s^p(j)
      newA(i,j) = pol(a, p(j));
     end; % if i
    end;
   end;

   % Newton recursion

   if verbose,
    disp('SPF: Start Newton iterative process.');
   end;

   stop = 0; it = 0; nwarn = 0;
   normB = norm(B);
   while ~stop,

    it = it + 1;
    A = newA;
    up = norm(A'*A-B)/normB; % relative residue

    % 1st stopping criterion : when B is factorized
    stop = (up < tolfact);

    if verbose,
     if stop,
      disp(['SPF: Iteration #' int2str(it) '  Residue =' num2str(up)]);
     else
      disp(['SPF: Iteration #' int2str(it) '  Residue =' num2str(up) ...
            ' (> ' num2str(tolfact) ')']);
     end;
    end;

    if ~stop,
     % solution with columnwise leading coefficient matrix triangular
     if strcmp(method, 'red'),
      X = axxab(A, 2*B, tol, method);
     else
      X = axxab(A, 2*B, tol, 'tri', method);
     end;
     if any(any(isnan(X))),
      disp('SPF: The symmetric polynomial equation has no solution.');
      error('The factorization failed !');
     end;

     newA = pzer((A+X)/2, tolzero);

     % A becomes unstable
     rootA = roots(newA); 
     maxroot = max(real(rootA));
     if maxroot > tolzero, 
      warning(['SPF: Unstable factor. Real part = ' num2str(maxroot) ' > 0.']);
      nwarn = nwarn + 1;
     end;
     
     if (nwarn > 10) | (it > 50),
      error('The iterative process does not seem to converge.');
     end;

    end;

   end;

   % computation of J : zero diagonal elements correspond to zero rows in A
   J = eye(n);
   for i = 1:n,
     if norm(A(i,:)) < tolzero,
       J(i,i) = 0;
     end;
   end; 

  else
   
   error('Invalid input argument.');

  end;  % if algorithm

  % multiplication by the unimodular matrix
  % corresponding to the diagonal reduction

  A = pzer(A*Udiagred, tolzero);
  pprop(A, typeB); 

  % Final check

  if verbose,
   if strcmp(method,'ext') | strcmp(method,'are')
    B0 = B{0}; A0 = A{0};
    Eps = J-(A0/B0)*A0';
    Eps = norm(Eps,1)/norm(B0,1);
   else
    Eps = B - A'*J*A;
    Eps = norm(Eps) / norm(B);
   end
   if isinf(Eps),
    disp('SPF: Final check is not performed.')
   elseif Eps > tolzero,
    warning(['SPF: Final check. High residue =' num2str(Eps)]);
   else
    disp('SPF: Final check is OK.');
   end;
  end;
 end;

else

 % #################
 % DISCRETE-TIME
 % #################

 if isempty(method), method = 'syl'; end;

 J = eye(n);

 if strcmp(type, 'forward'),
  B = B.';
 end;

 if issingular(B, tolzero),

  error('Matrix is singular.');

 elseif strcmp(method, 'syl') | strcmp(method, 'red'),

  if verbose,
   disp('SPF: Check positiveness of input matrix.');
  end;

  % initial matrix A0

  if verbose,
   disp('SPF: Build initial factor.');
  end;

  Bs0 = B{degB/2}; % absolute coefficient matrix

  if min(eig(Bs0)) < tolzero,
   error('The input matrix is not positive-definite on the unit circle.');
  end;

  newA = pol(chol(Bs0), 0);

  % Newton recursion

  if verbose,
   disp('SPF: Start Newton iterative process.');
  end;

  stop = 0; it = 0; nwarn = 0;
  normB = norm(B);
  while ~stop,

   it = it + 1;
   A = newA;
   up = norm(A'*A-B)/normB; % relative residue
   stop = (up < tolfact); % stop = 1 if the process converges

   if verbose,
    if stop,
     disp(['SPF: Iteration #' int2str(it) '  Residue =' num2str(up)]);
    else
     disp(['SPF: Iteration #' int2str(it) '  Residue =' num2str(up) ...
           ' (> ' num2str(tolfact) ')']);
    end;
   end;

   if ~stop,

    X = axxab(A, 2*B, tol, method, 'tri');
    newA = pzer((A+X)/2, tolzero);
    rootA = roots(newA);
    minroot = min(abs(rootA));
    if minroot < 1 - tolzero,
     warning(['SPF: Unstable factor. Module = ' num2str(minroot) ' < 1.']);
     nwarn = nwarn + 1;
    end;
 
    if (nwarn > 10) | (it > 50),
     error('The iterative process does not seem to converge.');
    end;

   end;

  end;

  % computation of J : zero diagonal elements correspond to zero rows in A
  J = eye(n);
  for i = 1:n,
   if norm(A(i,:)) < tolzero,
    J(i,i) = 0;
   end;
  end; 

  A = pzer(A, tolzero);

  if strcmp(type, 'forward'),
   B = B.';
   % column-wise transposition to avoid insertion
   % of zeros at the origin
   newA = pol(zeros(n));
   for i = 1:n,
    newA(:, i) = (A(:, i)').';
   end;
   A = newA;
  end;

  pprop(A, typeB);

  % Final check

  if verbose,
   residue = norm(B - A'*J*A);
   if residue > tolzero,
    warning(['SPF: Final check. High residue =' num2str(residue)]);
   else
    disp('SPF: Final check is OK.');
   end;
  end;

 else
 
  error('Invalid command option for discrete time case.');

 end;

end;

%end .. spf
