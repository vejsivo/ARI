function Z = null(A, varargin)
%NULL  Nullspace of polynomial
%
% The commmand  
%    Z = NULL(A) 
% computes a polynomial basis for the right null space of the
% polynomial matrix A, that is,
%    AZ = 0.
% If A has full column rank then Z is an empty polynomial matrix.
%
% If A is a pencil matrix, i.e. if the degree of A is equal to one,
% then the macro computes Z with a specialized algorithm taking advantage
% of the pencil structure of A.
%
% The commmand
%    Z = NULL(A,DEGREE) 
% seeks a basis of given degree DEGREE. If DEGREE is not specified
% then a miminal polynomial basis is computed by an iterative scheme. 
% The basis is minimal in the sense of Forney (see MINBASIS.) If 
% DEGREE is negative then the function directly computes an upper 
% bound degree basis.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: POL/RANK, XAB.

%   Author: D. Henrion, October 27, 1998.
%   Copyright 1998-2000 by Polyx, Ltd.
%   Modified by D. Henrion, February 2, 2000.
%               J. Jezek, July 19, 2001  arg checking
%               J. Jezek, Feb  28, 2003  arg checking

global PGLOBAL;

% verbose level
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');
           
% Parse optional input arguments.

tol = [];
degree = [];
if nargin > 1,
 for i = 1:length(varargin),
  arg = varargin{i};
  if ~isempty(arg),
   if isa(arg, 'double'), % matrix or scalar
    if ~any(size(arg) - [1 1]), % scalar
     if floor(arg) ~= arg, % non integer = tolerance
      if ~isreal(arg) | arg<0,  % here, exceptionally, arg>1
            % is not prohibited. It caused troubles. J.Jezek.
       error('Invalid tolerance.');
      else
       tol = arg;
      end;
     else % integer = degree
      degree = arg;
     end;
    else
     error(['Invalid ',nth(i+1),' argument.']);
    end;
   end;
  end;
 end;
end;

A = pol(A);  Avar = A.v; Ah = A.h;
if any(any(isnan(A))) | any(any(isinf(A))),
 error('Polynomial is not finite.');
end
rA = A.s(1); cA = A.s(2); degA = A.d;

% tolerance for rank computation
minA = min(abs(nonzeros(A.c)));
if isempty(tol),
 tolrank = PGLOBAL.ZEROING * minA;
else
 tolrank = tol * minA;
end;

% tolerance for zeroing
if isempty(tol),
 tolzero = PGLOBAL.ZEROING;
else
 tolzero = tol;
end;

% Null space extraction.

if isempty(A) | (degA < 0), % empty or zero polynomial matrix

 % null space is identity matrix

 Z = pol(eye(cA), 0);

elseif degA == 0, % constant polynomial matrix

 % standard null space through singular value decomposition

 if verbose,
  disp('NULL: Constant null space computation.');
 end;

 Z = pol(null(A.c), 0);

elseif degA == 1, % null space of a pencil matrix

 % The algorithm is from T.G.J. Beelen, G.W. Weltkamp "Numerical computation
 % of a coprime factorization of a transfer function matrix" Systems and
 % Control Letters, Vol. 9, pp. 281-288, 1987. It is based on the
 % Kronecker's canonical form of a singular form, obtained with the
 % algorithm described in P. Van Dooren "The computation of Kronecker's
 % canonical form of a singular pencil" Linear Algebra and its Applications,
 % Vol. 27, pp. 103-140, 1979.
 
 P = A;

 [n,m] = size(P);
 degP = deg(P);
 varP = pol([0 1],1,P.v);

 verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

 if isempty(tol),
  tol = PGLOBAL.ZEROING;
 end;

 A = polyval(P,0);
 B = polyval(P,1)-A;

 % Compute the null space of the pencil in two stages:
 % 1. Reduce the (singular) pencil s*B+A to generalized Schur form
 %    using Van Dooren's algorithm.
 % 2. Build the minimal degree polynomial basis for the null space
 %    using Beelen's algorithm.

 % First step: reduction of the singular pencil

 if verbose,
  disp('NULL: Reduction to Kronecker canonical form.');
 end;

 oldA = A; oldB = B;

 U = eye(n); V = eye(m); nn = max(m+1,n+1);
 row = zeros(1,nn+1); col = zeros(1,nn+1); step = 1;
 nrow = zeros(1,nn+1); ncol = zeros(1,nn+1);
 r = m; s2 = n;

 s = 1;
 while s > 0,

   % column compression
   nr = s2; nc = r;

   nrow(step) = nr; ncol(step) = nc;
   [B2,V2,r] = colcomp(B(1:s2, 1:r), tol); s = nc - r;
   row(step) = s; 

   if verbose,
    fprintf('NULL: Col. compr. (%dx%d) ', nr, nc);
   end;

   if s > 0,

    B(1:nr, 1:r) = B2(:, 1:r);
    B(1:nr, r+1:nc) = zeros(nr, s);
    B(nr+1:n, 1:nc) = B(nr+1:n, 1:nc) * V2;
    A(:, 1:nc) = A(:, 1:nc) * V2;
    V(:, 1:nc) = V(:, 1:nc) * V2;

    % row compression
    nr2 = nr; nc2 = s;
    [A2,U2,r2] = rowcomp(A(1:nr2, r+1:nc), tol); s2 = nr2 - r2;
    col(step) = r2;

    if verbose,
     fprintf(' Row compr. (%dx%d)', nr2, nc2);
    end;

    A(s2+1:nr2, r+1:nc) = A2(s2+1:nr2, :);
    A(1:s2, r+1:nc) = zeros(s2, nc2);
    A(1:nr2, 1:r) = U2 * A(1:nr2, 1:r);
    B(1:nr2, 1:r) = U2 * B(1:nr2, 1:r);
    U(1:nr2, :) = U2 * U(1:nr2, :);

   end;

   step = step + 1;

   if verbose,
    residue = norm(U*oldB*V-B) + norm(U*oldA*V-A);
    fprintf(' Residue = %f\n', residue);
   end;

  end;

  step = step - 2;

  % check if pencil has full column rank
  if step == 0,
  
   Z = pol(zeros(size(A,2),0));

  elseif (step > 0),

  % additional row index to ease block indexing
  if (col(step) > 0),
   step = step + 1;
   row(step) = ncol(step);
   cc = 0;
  end;
  rr = nrow(step+1); cc = ncol(step+1);

 U = flipud(U);
 V = fliplr(V);
 A = flipud(fliplr(A));
 B = flipud(fliplr(B));

 if verbose,
  residue = norm(U*oldB*V-B) + norm(U*oldA*V-A);
  fprintf('NULL: Final residue = %f\n', residue);
 end;

 % Second step: computation of the null space basis
  
 N = cell(step,step); AI = cell(1,step); dim = zeros(1,step);
 rowind = cumsum([1 row]); colind = cumsum([1 col]);

 if verbose,
  fprintf('NULL: Null spaces and right inverses of block diagonal terms.\nNULL: Block ');
 end;

 for j = 1:step,

  if verbose,
   fprintf('%d ', j);
  end;

  AB = A(colind(j):colind(j+1)-1, rowind(j):rowind(j+1)-1);

  [N{j,j} AI{j} dim(j)] = right(AB, tol);

  if dim(j) ~= row(j)-col(j),
   error('NULL: Inconsistent null space dimensions.');
  end;

 end;

 if verbose,
  fprintf('\nNULL: Iterative computation of off-diagonal terms.\n');
  fprintf('NULL: Column ');
 end;
 
 for j = 2:step,

  if dim(j) > 0,

  if verbose,
   fprintf('%d ', j);
  end;

   for i = j-1:-1:1,

    BA = B(colind(i):colind(i+1)-1, rowind(i+1):rowind(i+2)-1);
    AA = A(colind(i):colind(i+1)-1, rowind(i+1):rowind(i+2)-1);
    M = (varP*BA+AA)*N{i+1,j};
    for k = i+2:j,
     BA = B(colind(i):colind(i+1)-1, rowind(k):rowind(k+1)-1);
     AA = A(colind(i):colind(i+1)-1, rowind(k):rowind(k+1)-1);
     M = M + (varP*BA+AA)*N{k,j};
    end;

    N{i,j} = AI{i}*M;

   end;

  end;

 end;

 % Build minimal polynomial basis

 if verbose,
  fprintf('\n');
  disp('NULL: Build minimal polynomial basis.');
 end;

 cZ = sum(dim); rZ = sum(row)+cc;
 dimind = cumsum([1 dim]);

 Z = pol(zeros(rZ, cZ));
 for i = 1:step,
  for j = i:step,
   Z = assign(Z, rowind(i):rowind(i+1)-1, dimind(j):dimind(j+1)-1, N{i,j});
  end;
 end;

 Z = pzer(V*Z);

 if verbose,

  % Final check
  residue = norm(P*Z);
  disp(['NULL: Final check. Residue = ' num2str(residue)]);

 end;

 end;
 
 if verbose,
  if isempty(Z),
   disp('NULL: No polynomial right null space');
  end;
 end;

else, % polynomial null space

 % 1. Evaluate the rank of polynomial matrix A with overloaded macro RANK.
 % 2. With CEF and NULLREF, extract standard null spaces of Sylvester matrices
 %    of A of increasing orders. Stop when the number of columns of the
 %    polynomial null space is equal to the nullity of A.

 rkA = rank(A, tolrank);
 
 nullA = cA - rkA; % nullity of A

 if nullA > 0,

  % upper bound on degree of null space, cf. [1]
  cdeg = deg(A, 'col'); cdeg(cdeg < 0) = 0; cdeg = -sort(-cdeg);
  rdeg = deg(A, 'row'); rdeg(rdeg < 0) = 0; rdeg = -sort(-rdeg);
  rdeg = [rdeg; zeros(cA-rA, 1)];
  degZmax = min(sum(cdeg(1:(cA-1))), sum(rdeg(1:(cA-1))));

  if ~isempty(degree) & (degree > degZmax),
   if verbose,
    disp(['NULL: User-supplied degree greater than ' ...
          'maximum expected degree (' int2str(degZmax) ').']);
   end;
  end;  

  % CEF and NULLREF work row-wise, so matrix A must be transposed
  A = A.';

  % order of Sylvester matrix
  if isempty(degree),
   degZ = 0;
   if verbose,
    disp('NULL: Seek minimum degree polynomial null space.');
   end;
  elseif degree >= 0,
   degZ = degree; degZmax = degree;
   if verbose,
    disp(['NULL: Seek polynomial null space of degree ' int2str(degZ) '.']);
   end;
  else,
   degZ = degZmax;
   if verbose,
    disp(['NULL: Seek polynomial null space of degree ' int2str(degZ) '.']);
   end;
  end;

  % Main loop: computation of null spaces of
  % Sylvester matrices of increasing order.

  delta = degZmax / 10; step = 1;
  degZ = degZ - 1;
  ndep = 0;

  while (degZ < degZmax) & (ndep < nullA),

   % Order of Sylvester matrix: the step is increased
   % to accelerate the reduction.
   degZ = min(degZmax, degZ + floor(step));
   step = step + delta;

   if verbose,
    fprintf('NULL: Degree %3d (Max %3d) ', degZ, degZmax);
   end;

   % Transform Sylvester matrix into column echelon form.

   RA = sylv(A, degZ);
   [Z, irow] = cef(RA, tol);

   % Number of dependent rows.

   % dependent row indices in Sylvester matrix
   % transformed into dependent row indices in polynomial matrix
   drow = 1:size(RA,1); drow(irow) = 0; drow = rem(drow-1,cA)+1;
   % number of dependent rows in polynomial matrix
   count = zeros(1, cA); count(drow(drow > 0)) = 1;
   ndep = sum(count);

   if verbose,
    disp(['  Number of missing rows = ' int2str(nullA - ndep)]);
   end;

  end; % while
 
  if nullA <= ndep,

   % Compute coefficients of null space.
   Z = nullref(Z, irow);

   % Extract primary dependent rows.

   dep = zeros(nullA, 1); j = 1;
   drow = drow(drow > 0);
   for i = 1:cA,
    if count(i),
     ind = find(drow == i);
     dep(j) = ind(1); j = j + 1;
    end;
   end;
   dep = dep(1:nullA);

   % Scale each row
   for i = dep',
    Z(i, :) = Z(i, :) / norm(Z(i, :));
   end;
   Z = pol(Z(dep, :), degZ).';

   % zeroing
   Z = pzer(Z, tolzero);

   if verbose,
    disp(['NULL: A polynomial null space of degree ' ...
          int2str(deg(Z)) ' was found.']);
   end;

  else % polynomial null space of given degree does not exist

   Z = pol(zeros(cA,0));
  
   if verbose,
    disp(['NULL: No polynomial null space of degree ' ...
          int2str(degZ) ' was found.']);
   end;

  end;

 else % full column rank = no polynomial null space

  Z = pol(zeros(cA,0));

  if verbose,
   disp('NULL: No polynomial null space.');
  end;

 end;

end;
props(Z,Avar,Ah);

% SUBROUTINES USED FOR THE PENCIL ALGORITHM

% Column compression
% A*V = B = [B1 0], r = size(B1,2)
function [B,V,r] = colcomp(A, tol)
if ~isempty(A),
 [U,S,V] = svd(A);
 B = U*S;
 if min(size(S)) > 1, S = diag(S); end;
 r = sum(S > max(size(A))*tol);
else
 B = zeros(size(A)); V = eye(size(A,2)); r = size(A,2);
end;

% Row compression
% U*A = B = [0; B2], r = size(B2,1)
function [B,U,r] = rowcomp(A, tol)
if ~isempty(A),
 [U,S,V] = svd(A);
 U = U'; B = S*V';
 if min(size(S)) > 1, S = diag(S); end;
 r = sum(S > max(size(A))*tol);
 B = B([(r+1):size(A,1) 1:r],:);
 U = U([(r+1):size(A,1) 1:r],:);
else
 B = zeros(size(A)); U = eye(size(A,1)); r = size(A,1);
end;

% Right null space and right inverse
% A*Z = 0, A*R = -1, n = size(Z,2)
function [Z,R,n] = right(A, tol)
if ~isempty(A),
 [U,S,V] = svd(A);
 if min(size(S)) > 1, S = diag(S); end;
 r = sum(S > max(size(A))*tol);
 S = diag(ones(r,1)./S(1:r));
 Z = V(:,r+1:size(V,2));
 n = size(Z,2);
 R = -V(:,1:r)*S*U(:,1:r)';
else
 Z = eye(size(A,2)); R = A'; n = size(Z,2);
end;

% Polynomial matrix assignment without using subscript references
% A(rowind,colind) = B
function A = assign(A,rowind,colind,B)
if ~isempty(B),
 A = pol(A); B = pol(B);
 Ac = A.c; Bc = B.c;
 [a1,a2,a3] = size(Ac); [b1,b2,b3] = size(Bc);
 if ~a3, Ac = zeros(a1, a2); a3 = 1; end;
 if ~b3, Bc = zeros(b1, b2); b3 = 1; end;
 Ac = cat(3, Ac, zeros(a1, a2, b3-a3));
 Bc = cat(3, Bc, zeros(b1, b2, a3-b3));
 Ac(rowind, colind, :) = Bc;
 A = pol(Ac(:, :), max(a3, b3) - 1);
end;

%end .. @pol/null
