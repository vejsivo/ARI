function [T,U,rowind] = tri(A, varargin)
%TRI  Triangular or staircase form of a polynomial matrix
%
% Given an arbitrary polynomial matrix A, the commands
%    [T,U,ROWIND] = TRI(A)
%    [T,U,ROWIND] = TRI(A,'col')
% compute a column staircase form of A. In particular, if A is 
% nonsingular then T is a lower-left triangular matrix. The unimodular 
% reduction matrix U is such that
%     AU = T
% The row vector ROWIND is such that ROWIND(I) is the row index of
% the uppermost nonzero entry in the I-th column of T. Note that 
% ROWIND(I) = 0 if there are no nonzero entries in I-th column of T.
% As a result, the number of non-zero elements in ROWIND is the
% algorithm's idea of the rank of A.
%
% Off-diagonal elements in the staircase form are not necessarily reduced,
% so the output generally differs from that of the macro HERMITE.
%
% Similarly, 
%   [T,U,COLIND] = TRI(A,'row') 
% computes a row staircase form of A, the unimodular reduction matrix U 
% such that UA = T, and a column vector COLIND such that COLIND(I) is the 
% column index of the leftmost non-zero entry in the I-th row of T.
%
% The commmand
%    TRI(A,DEGREE) 
% seeks a reduction matrix of given degree DEGREE. If DEGREE is not 
% specified then a reduction matrix of minimum overall degree is computed 
% by an iterative scheme. If DEGREE is negative then a reduction matrix 
% of maximum achievable degree is returned.
%
% The commmand
%    TRI(A,'syl') 
% performs the transformation through stable reductions of Sylvester matrices. 
% This method is preferable numerically. It is the default method.
%
% The commmand
%    TRI(A,'gau') 
% performs the transformation through a modified version of Gaussian elimination. 
% This method is preferable esthetically.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: HERMITE.

%    Author: D. Henrion,  November 16, 1998.
%    Updated to 3.0 by D. Henrion, August 30, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.
%    $ Revision $   $ Date 17-Jul-2001   J.Jezek  $
%                   $ arg checking, sampling period, empty case  $

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin<1,
   error('Not enough input argments.');
end;
eval('A = pol(A);', 'error(peel(lasterr));');

if any(any(isnan(A))) | any(any(isinf(A))),
 error('Polynomial is not finite.')
end;

% verbose flag
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% Default options.

tol = [];
method = 'syl';
form = 'col';
degree = [];

% Parse optional input arguments.

invalid = 0;
if nargin > 1,
 for i = 1:length(varargin),
  arg = varargin{i};
  if ~isempty(arg),
   if isa(arg, 'char'), % valid strings
    if strcmp(arg, 'syl'),
     method = 'syl';
    elseif strcmp(arg, 'gau'),
     method = 'gau';
    elseif strcmp(arg, 'col'),
     form = 'col';
    elseif strcmp(arg, 'row'),
     form = 'row';
    else,
     error('Invalid command option.');
    end;     
   elseif isa(arg, 'double'), % matrix or scalar
    if ~any(size(arg) - [1 1]), % scalar
     if floor(arg) ~= arg, % non integer = tolerance
      tol = arg;
     else % integer = degree
      degree = arg;
     end;
    else, % matrix argument
     error(['Invalid ',nth(i+1),' argument.']);
    end;
   else,
    error(['Invalid ',nth(i+1),' argument.']);
   end;
  end;
 end;
end;

if ~isempty(tol),
   if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;   
if strcmp(method, 'gau') & ~isempty(degree),
 error(['The degree of the reduction matrix cannot be specified ' ...
       'for Gaussian elimination.']);
end;

if strcmp(form, 'row'),
 A = A.';
end;

[rA,cA] = size(A); dA = deg(A);
var = A.var; h = A.h;

% tolerance for zeroing
me = 1;
if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else 
 tolzero = tol * me;
end;

if isempty(A) | dA < 0,

 % Empty matrix  or  zero matrix.
 T = pol(zeros(rA, cA));
 U = pol(eye(cA));
 rowind = zeros(1, cA);

elseif strcmp(method, 'syl'),

 % ********************************
 % SYLVESTER MATRIX METHOD, cf. [1]
 % ********************************

 if dA == 0,
 
  % triangularization of a constant matrix
  %---------------------------------------

  % If the reduction matrix is required,
  % append an identity matrix at the bottom of A.

  if nargout > 1,
   A = [A.coef; eye(cA)];
  else
   A = A.coef;
  end;

  % Call CEF, the constant matrix version of the present macro. 

  [T,rowind] = cef(A, tol);

  rowind(rowind > rA) = 0;
  rowind = [rowind zeros(1,cA-length(rowind))]; 

  % Recover reduction matrix.

  if nargout > 1,
   U = T(rA+1:rA+cA, :);
   U = pzer(pol(U), tolzero);
  end; 

  T = T(1:rA, :);
  T = pzer(pol(T), tolzero);

 else,

  % triangularization of a polynomial matrix
  %-----------------------------------------

  % upper bound on the degree of the reduction matrix
  cdeg = deg(A, 'col'); cdeg(cdeg < 0) = 0; cdeg = -sort(-cdeg);
  rdeg = deg(A, 'row'); rdeg(rdeg < 0) = 0; rdeg = -sort(-rdeg);
  rdeg = [rdeg; zeros(cA-rA, 1)];
  maxdegree = min(sum(cdeg(1:(cA-1))), sum(rdeg(1:(cA-1))));

  if ~isempty(degree) & (degree > maxdegree),
   maxdegree = degree;
   warning(['User-supplied degree greater than ' ...
          'maximum expected degree (' int2str(maxdegree) ').']);
  end;  
 
  if verbose,
   disp('TRI: Sylvester matrix method.');
   if ~isempty(degree),
    if degree >= 0,
     disp(['TRI: Seek reduction matrix of degree ' int2str(degree) '.']);
    else
     disp('TRI: Seek maximum degree reduction matrix.');
    end;
   else
    disp('TRI: Seek minimum degree reduction matrix.');
   end;
  end;

  % Append identity matrix to recover reduction matrix.

  A = [A; eye(cA)];
  rA = rA + cA;

  % Permute matrix coefficients of A.

  Ac = A.coef;
  r = rA*(dA+1);
  B = zeros(r , cA);
  for i = 0:dA,
   B(1+i*rA : (i+1)*rA, :) = Ac(:, 1+(dA-i)*cA : (dA-i+1)*cA);
  end;

  if isempty(degree),
   dU = 0;
  elseif degree < 0,
   dU = maxdegree;
  else
   dU = degree;
  end;

  lowerbound = dU;
  upperbound = maxdegree;

  solution = 0; lost = 0; givendegree = ~isempty(degree); stay = 1;

  iteration = 0;

  % Binary search for minimum degree reduction matrix.
  % --------------------------------------------------

  while stay,

   if iteration,
    dU = ceil((upperbound + lowerbound) / 2);
   end;

   iteration = iteration + 1;

   if verbose,
    fprintf('TRI: Attempt #%3d  Degree %3d (Min %3d / Max %3d) ', ...
     iteration, dU, lowerbound, upperbound);
   end;

   % build row permuted Sylvester matrix
   %---------------------------------------

   d = dA+dU+1;
   cmax = cA*(dU+1); rmax = rA*d;
   Res = zeros(rmax, cmax);

   index1 = zeros(1, r); index2 = zeros(1, r);
   k = 0; l = 0; for i = 1 : rA, for j = 0 : dA,
    k = k+1; l = l+1;
    index1(k) = i + j*rA; index2(k) = l;
   end; l = l+dU; end;

   for i = 0 : dU,
    Res(i+index2, 1+cA*i : cA*(i+1)) = B(index1, :);
   end;

   % triangularize permuted Sylvester matrix
   %---------------------------------------------------

   [T,irow] = cef(Res, tol);

   % retrieve triangular shape

   if length(irow) < cmax,
    if verbose, fprintf('*\n'); end;
    disp('Sylvester matrix is rank deficient. Triangularization failed.');
    error('Try to modify the tolerance.');
   end;

   shape = zeros(1, rA);
   for i = 1:cmax,
    shape(floor((irow(i)-1)/d)+1) = i; 
   end;

   % count leading entries
   nentries = sum(shape & shape);

   triangular = 0;
   if nentries > cA,
    if verbose, fprintf('*\n'); end;
    disp('Sylvester matrix is not reduced. Triangularization failed.');
    error('Try to modify the tolerance.');
   elseif nentries == cA,
    triangular = 1;
    solution = 1;
   end;

   % Next degree.
   % ------------
 
   if triangular,

    upperbound = dU; 
    oldRes = Res; oldT = T; oldshape = shape;
    lost = 0; 

    if verbose,
     fprintf('Yes.\n');
    end;

   else

    lowerbound = dU; lost = 1;

    if verbose,
     fprintf('No.\n');
    end;

   end;

   if givendegree, givendegree = 0; stay = 0; end; % given degree: only once
   if ~solution & (dU == maxdegree - 1), stay = 1; % last chance
   else stay = (upperbound - lowerbound > 1); end;

  end; % binary search

  if ~solution,
   disp('Triangularization failed.');
   error('Try to modify the tolerance.');
  end;

  if lost,
   dU = upperbound; Res = oldRes; T = oldT; shape = oldshape;
  end;

  append = 0;

  % Recover triangular form.
  % ------------------------ 

  d = dA+dU+1; rA = rA - cA;
  index = zeros(1, rmax);
  k = 0; for i = 1 : d, for j = 0 : rA+cA-1,
   k = k+1; index(k) = i + j*d;
  end; end;

  S = T(index, shape(shape > 0));

  T = zeros(rA, cA*d);
  for i = 1:d,
   T(:, 1+(i-1)*cA : i*cA) = S(1+(d-i)*(rA+cA) : rA+(d-i)*(rA+cA), :);
  end;
  T = pzer(pol(T, d-1, var), tolzero);
  T.h = h;

  if verbose,
   disp(['TRI: Staircase form of degree ' int2str(deg(T)) '.']);
  end;

  if nargout > 1,

   % Recover reduction matrix.
   % -------------------------

   U = zeros(cA, cA*(dU+1));
   for i = 1:dU+1,
    U(:, 1+(i-1)*cA : i*cA) = S(1+rA+(d-i)*(rA+cA) : rA+cA+(d-i)*(rA+cA), :);
   end;
   U = pzer(pol(U, dU, var), tolzero);
   U.h = h;

   if verbose,
    disp(['TRI: Reduction matrix of degree ' int2str(deg(U)) '.']);
   end;

  end;

  % Row indices of leading entries.
  % -------------------------------

  rowind = (shape & shape) .* (1:(rA+cA));
  rowind = rowind(rowind > 0); rowind(rowind > rA) = 0;

 end; % if constant

else

 % ********************
 % GAUSSIAN ELIMINATION
 % ********************

 if verbose,
  disp('TRI: Gaussian elimination method.');
  disp('TRI: Seek minimum degree reduction matrix.');
 end;

 [rA, cA] = size(A); dA = deg(A);

 % upper bound on the degree of the reduction matrix
 cdeg = deg(A, 'col'); cdeg(cdeg < 0) = 0; cdeg = -sort(-cdeg);
 rdeg = deg(A, 'row'); rdeg(rdeg < 0) = 0; rdeg = -sort(-rdeg);
 rdeg = [rdeg; zeros(cA-rA, 1)];
 maxdegree = min(sum(cdeg(1:(cA-1))), sum(rdeg(1:(cA-1))));

 % Append identity matrix to A.

 A = [A; eye(cA)]; rA = rA+cA;

 % Permute matrix coefficients of A
 % to build row permuted Sylvester matrix.

 r = rA*(dA+1);
 Ac = A.coef;
 A = zeros(r , cA);
 index = zeros(1, r);
 k = 0; for i = 1 : rA, for j = 0 : dA,
  k = k+1; index(k) = i + j*rA;
 end; end;
 for i = 0:dA,
  A(1+i*rA : (i+1)*rA, :) = Ac(:, 1+(dA-i)*cA : (dA-i+1)*cA);
 end;
 A = A(index, :);

 pivot = zeros(1, cA);
 found = zeros(1, rA);
 entry = zeros(r, rA);

 dU = 0;

 % Loop until all the leading entries are found
 % or maximum reduction degree matrix is reached
 % --------------------------------------------

 while sum(found) < cA,

  if verbose,
   fprintf(['TRI: Degree = %3d  ' ...
            'Number of missing leading entries = %3d\n'], dU, cA-sum(found));
  end;

  rmax = (dA+dU+1)*rA;

  r = 1; c = 1;

  if dU > 0,
  
   % build shifted permuted Sylvester matrix

   index = zeros(1, rmax);

   k = 0; l = 0; for i = 1 : rA, for j = 1 : dA+dU,
    k = k+1; l = l+1; index(k) = l;
   end; l = l+1; end;
   for i = 1 : rA,
    index(rmax - rA + i) = i*(dA+dU+1);
   end;

   B1 = [A(:, 1:dU*cA); zeros(rA, dU*cA)];
   B2 = [A(:,1+(dU-1)*cA:dU*cA); zeros(rA, cA)];

   A(index, 1:dU*cA) = B1;
   index = index + 1; index(rmax) = 1;
   A(index, 1+dU*cA:(dU+1)*cA) = B2;

   entry(index, :) = [entry; zeros(rA)];

   % update pivots

   for i = 1:dU*cA,
    pivot(i) = pivot(i) + floor((pivot(i)-1)/(dA+dU));
   end;
   pivot = [pivot zeros(1, cA)];

   c = dU*cA + 1;

  end;

  % Inspect shape
  % -------------

  me = max(size(A)) * norm(A, 'inf');
  if isempty(tol),
   toltri = PGLOBAL.ZEROING * me;
  else
   toltri = tol * me;
  end;

  while c <= (dU+1)*cA,

   % seek non-zero entry in column
   next = 1;
   while next,
    if r <= rmax,
     if norm(A(r, c)) < toltri, r = r+1;
     else next = 0; end;
    else next = 0; end;
   end;

   % is it a new pivot ?

   if r > rmax,

     disp('Some leading entries are missing.');
     error('The triangularization failed. Try to modify the tolerance.');

   elseif any(pivot == r),

     % no there is already a pivot in the same row
     % cancel column
     colpivot = find(pivot == r);
     A(r:rmax, c) = A(r:rmax, colpivot) - A(r:rmax, c) / A(r, c);
     A(r, c) = 0;

     % look for a pivot in the same column
     if r < rmax,
      r = r + 1;
     end;

   else
   
     % yes it is a pivot, so store it
     pivot(c) = r;  

     % normalize column
     A(r+1:rmax, c) = A(r+1:rmax, c) / A(r, c);
     A(r, c) = 1;

     % store column
     row = 1+floor((r-1)/(dA+dU+1));
     found(row) = 1;
     entry(:, row) = A(:, c);
    
     % next column, top row
     c = c+1; r = 1;

   end;

  end; % while c <= (dU+1)*cA

  if dU > maxdegree,
   disp('Maximum reduction degree was reached and some leading entries are missing.');
   error('The triangularization failed. Try to modify the tolerance.');
  else
   dU = dU+1;
  end;

 end; % until all leading entries are found

 % Recover unimodular reduction matrix U
 % -------------------------------------

 entry = entry(:, find(found));
 dU = dU - 1; rA = rA - cA;
 T2 = entry(1 : rA*(dA+dU+1), :);
 U2 = entry(1+rA*(dA+dU+1) : rmax, :);

 if nargout > 1,

  index = zeros(1, cA*(dA+dU+1));
  k = 0; for i = 1 : cA, for j = 0 : dA+dU,
    k = k+1; index(k) = i + j*cA;
  end; end;
  U2(index, :) = U2;
  U = zeros(cA, cA*(dU+1));
  for i = 0 : dU,
   U(:, 1+(dU-i)*cA:(dU-i+1)*cA) = U2(1+(i+dA)*cA:(i+dA+1)*cA, :);
  end;

  U = pzer(pol(U, dU, var), tolzero);
  U.h = h;

  if verbose,
   disp(['TRI: Reduction matrix of degree ' int2str(deg(U)) '.']);
  end;

 end;

 % Recover triangular form T.
 % --------------------------

 index = zeros(1, rA*(dA+dU+1));
 k = 0; for i = 1 : rA, for j = 0 : dA+dU,
   k = k+1; index(k) = i + j*rA;
 end; end;
 T2(index, :) = T2;
 T = zeros(rA, cA*(dA+dU+1));
 for i = 0 : dA+dU,
  T(:, 1+(dA+dU-i)*cA:(dA+dU-i+1)*cA) = T2(1+i*rA:(i+1)*rA, :);
 end;

 T = pzer(pol(T, dA+dU, var), tolzero);
 T.h = h;

 if verbose,
  disp(['TRI: Staircase form of degree ' int2str(deg(T)) '.']);
 end;

 % Row indices of leading entries.
 % -------------------------------

 rowind = found .* (1:(rA+cA));
 rowind = rowind(rowind > 0); rowind(rowind > rA) = 0;

end;

% ************
% LAST TUNINGS
% ************

if strcmp(form, 'row'),
 T = T.';
 if nargout > 1, U = U.'; end;
 if nargout > 2, rowind = rowind'; end;
end;

%end .. tri
