function AR = fact(A, Z, tol)
%FACT    Polynomial matrix factor extraction
%
% Given a non-singular polynomial matrix A and a vector Z containing 
% a subset of the zeros of A, the instruction
%        AR = FACT(A,Z)
% extracts from A a minimum degree, column and row reduced right
% polynomial matrix factor AR such that Z contains all the zeros of AR and
%        A = AL * AR.
% The left factor AL may be retrieved with the instruction AL = A/AR.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%    Author: D. Henrion, September 22, 1998.
%    Updated to 3.0 by D. Henrion, August 30, 2000.
%    Last modified by D. Henrion, November 20, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

%    The macro proceeds by interpolation. Right factor AR is extracted from
%    the null-space of a matrix built with the characteristic vectors of A
%    corresponding to the zeros in Z. Macro CHARACT is used for char. vector
%    computation. Left factor AL is retrieved with macro XAB.

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin<2,
   error('Not enough input arguments.');
end;
eval('A = pol(A);', 'error(peel(lasterr));');

[n, cA] = size(A);
if n ~= cA,
 error('Matrix is not square.');
end;

if ~isa(Z, 'double') | min(size(Z)) > 1,
   error('Invalid 2nd argument.');
end;

if nargin < 3 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

if rank(A, tol) < n,
 error('Matrix is singular.');
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% Number of zeros, including multiplicities

nzeros = length(Z);

% Detect multiple zeros

[Y,i] = sort(Z); Y = Z(i); Z = [];
m = length(Y); i = 1;
while i <= m,
 y = Y(i); j = 1;
 while (i + j <= m) & (abs(y - Y(i + j)) < tol), j = j + 1; end;
 Z = [Z; [y j]]; % zero and multiplicity
 i = i + j;
end;

% Number of distinct zeros

ndzeros = size(Z, 1);

% Compute characteristic vectors
% and interpolated matrix

if verbose,
 disp('FACT: Compute characteristic vectors and interpolated matrix.');
end;

S = mono([0:nzeros]'); % polynomial basis
Interpol = zeros(n*(nzeros+1), nzeros); % interpolated matrix
col = 1; % column in interpolated matrix

% For each distinct zero, build columns of interpolated matrix

for i = 1:ndzeros,

 zero = Z(i, 1);
 mult = Z(i, 2);

 % Compute characteristic vectors

 [vect, multmax, geom] = charact(A, zero, tol);

 if isempty(vect),

  warning([num2str(zero) ' is not a zero of the input polynomial matrix.']);

 else

  if mult > multmax,
   warning(['Multiplicity of ' num2str(zero) ' must not exceed ' ...
    int2str(multmax) '.']);
   mult = multmax;
  end;

  if verbose,
   disp(['FACT: Zero = ' num2str(zero) '  Mult. = ' int2str(mult) ...
        ' (Alg. = ' int2str(multmax) '  Geom. = ' int2str(geom) ').']);
  end;

  % Charac. vect. in interpolated matrix

  chain = 0; % char. vector chain index
  order = 1; % derivative order in chain

  while mult > 0,

   chain = chain + 1;
   if chain > geom,
    chain = 1; order = order + 1;
   elseif isempty(vect{chain, order}),
    chain = 1; order = order + 1;
   end;

   factorial = 1;
   for j = 0:order-1,
    Interpol(:, col) = Interpol(:, col) + ...
      kron(polyval(polyder(S, j), zero) / factorial, vect{chain, order-j});
    factorial = factorial * (j+1);
   end;

   col = col + 1; mult = mult - 1;

  end;

 end;

end;

% Extract right null-space of the interpolated matrix

if verbose,
 disp('FACT: Compute null-space of interpolated matrix.');
end;

[K, irow] = cef(Interpol, tol);
K = nullref(K, irow);
rK = size(K, 1);

% Indices of dependent rows

drow = 1:n*(nzeros+1); drow(irow) = 0; drow = drow(drow > 0);

% Compute column degrees:
% they correspond to column indices of leading entries in K

if verbose,
 disp('FACT: Compute column degrees in right factor.');
end;

rowsel = [];
coldegAR = zeros(1, n);

for i = 1:n,
 row = find(rem(drow-1, n)+1 == i);
 if length(row) < 0,
  error('FACT: Invalid column degrees. Try to modify the tolerance.');
 end;
 rowsel = [rowsel row(1)];
 coldegAR(i) = (drow(row(1)) - i) / n;
end;

if verbose,
 fprintf('FACT: Expected column degrees =');
 for i = 1:n,
  fprintf(' %2d', coldegAR(i));
 end;
 fprintf('.\n');
end;

degAR = max(coldegAR);

colsel = [];
for j = 0:nzeros, for i = 1:n,
 if coldegAR(i) >= j, colsel = [colsel i + j*n]; end;
end; end;

% Retrieve right factor AR

if verbose,
 disp('FACT: Right factor extraction.');
end;

coefAR = zeros(cA, cA*(degAR+1));
coefAR(:, colsel) = K(rowsel, colsel);
AR = pzer(pol(coefAR, degAR, A.var), tol);

%end .. fact

