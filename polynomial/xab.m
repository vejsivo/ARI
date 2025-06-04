function [X,K] = xab(A,B,varargin)
%XAB  Matrix polynomial equation solver
%
% The commmand
%    X0 = XAB(A,B) 
% finds a particular solution of the linear matrix polynomial equation
%    XA = B
% If no polynomial solution exists then all the entries in X0 are set 
% equal to NaN.
%
% The commmand
%    X0 = XAB(A,B,DEGREE) 
% seeks a solution X0 of degree DEGREE. If DEGREE is not specified then 
% a solution of minimum overall degree is computed by an iterative scheme.
% If DEGREE is negative then the function directly computes an upper bound 
% degree solution.
%
% The commmand
%    X0 = XAB(A,B,'sqz') 
% seeks a solution X0 with 'squeezed' column degrees. If N is the nullity 
% of A then the N rightmost column degrees are minimized, at the expense 
% of increasing the degrees in the other columns.
%
% If DEGREES is a vector of zeros and ones such that DEGREES(i) = 1 and
% DEGREES(j) = 0 then the function 
%    X0 = XAB(A,B,'sqz',DEGREES) 
% attempts to minimize the degree of the i-th column in X0, provided that 
% the degree of the j-th column may increase.
%
% The commmand
%    XAB(A,B,'syl') 
% solves the equation through the Sylvester matrix method. This is the 
% default method.
%
% The commmand
%    [X0, K] = XAB(A,B) 
% also computes the left null-space of A so that all the solutions 
% to XA = B may be parametrized as
%       X = X0 + TK
% where T is an arbitrary polynomial matrix.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: POL/NULL, AXB.

%    Author: D. Henrion, October 14, 1998.
%    Modified by J.Jezek, May 24, 2000.
%    Updated to 3.0 by D. Henrion, September 18, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.FORMAT;','painit;');

if nargin < 2,
 error('Not enough input arguments.');
end;

eval('A = pol(A); B = pol(B);', ...
   'error(peel(lasterr));');

if any(any(isnan(A))) | any(any(isnan(B))) | ...
   any(any(isinf(A))) | any(any(isinf(B))),
 error('Polynomial is not finite.');
end;

% Data type
[tv,typeA,A,B] = testvp(A,B);
if tv==2,
   error('Inconsistent variables.');
elseif tv==0,
   warning('Inconsistent variables.');
end;

[th,Xh,A,B] = testhp(A,B,typeA);
if th==0,
   warning('Inconsistent sampling periods.');
end;

% verbose level
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

[rA, cA] = size(A); degA = A.d;
[rB, cB] = size(B); degB = B.d;

scalarA = (rA == 1) & (cA == 1);
scalarB = (rB == 1) & (cB == 1);

if ~scalarA & (cA ~= cB),
 error('Matrices of inconsistent dimensions.');
end;

% Default options.

tol = [];
method = 'syl';
degree = [];
columns = []; % degree squeezing

% Handle the case when XAB is called from other equation solvers.

if nargin > 2,
 if isa(varargin{1}, 'cell'),
  varargin = varargin{1};
 end;
else
 varargin = [];
end;

% Parse optional input arguments.

invalid = 0;
lv = length(varargin);
if lv>0,
 for i = 1:lv,
  arg = varargin{i};
  if ~isempty(arg),
   if isa(arg, 'char'), % valid strings
    if strcmp(arg, 'syl'),
     method = 'syl';
    elseif strcmp(arg, 'sqz'),
     columns = zeros(1, rA);
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
    else, % matrix argument: vector of column indices
     if any(sort(size(arg)) - [1 rA]),
      error(['Invalid length of the degree index vector; must be ' int2str(rA) '.']);
     else,
      columns = arg;
     end;
    end;
   else,
    error(['Invalid ',nth(i+2),' argument.']);
   end;
  end;
 end;
end;

% tolerance for zeroing:
% relative to optional input parameter TOL and to elements in A and B

maxA = max(abs(nonzeros(A.c)));
minB = 1; minanBc = min(abs(nonzeros(B.c)));
if ~isempty(minanBc),
   minB = max(minB, minanBc);
end;
if maxA > 0, me = min(minB, 1/maxA);
else me = minB;
end;

if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else
 tolzero = tol * me;
end;

% relative tolerance for rank computation:

me = 1e-2; % one percent

if isempty(tol),
 tolrank = PGLOBAL.ZEROING * me;
else
 tolrank = tol * me;
end;

% *********************
% Handle special cases.
% *********************

solution = 1;

if isempty(A) | isempty(B),

 X = pol(zeros(rB,rA));
 solution = 1;
 
elseif (degA >= 0) & (degB < 0), % B is zero but not A

 if scalarA,
  verbose = 0;
  X = pol(zeros(rB, cB));
 else
  X = pol(zeros(rB, rA));
 end;

elseif scalarA & (degA <= 0), % A is a constant scalar

 verbose = 0;
 if degA < 0, % A is zero

  X = lcoef(B, 'ent');
  X = pol(Inf * X); X.h = Xh;

 else, % A is non-zero

  X = B.c / A.c;
  X = pol(X(:,:), degB, typeA); X.h = Xh;
  X = pzer(X, tolzero);
 
 end;

else

 if scalarA,
  if ~scalarB, % A is a non-zero polynomial scalar and B a matrix
   A = A * eye(cB);
   rA = cB; cA = cB;
  end;
 end;

 if strcmp(method, 'syl'),

  % -----------------------
  % SYLVESTER MATRIX METHOD
  % -----------------------

  % Necessary condition of existence of a solution.

  rkA = rank(A, tolrank);

  if rkA < cA,
   rkAB = rank([A; B], tolrank);
   if rkA < rkAB, % no solution
    solution = 0;
   end;
  end;

  if solution  % a solution may exist

   % Bounds on degree of X.
   % ......................

   if (degA > 0) | (degB > 0),
  
    % Lower bound.

    cdA = deg(A, 'col'); cdA(cdA < 0) = 0;
    cdB = deg(B, 'col'); cdB(cdB < 0) = 0;
    mindegree = max([cdB-cdA 0]); 

    if ~isempty(degree) & (degree >= 0) & (degree < mindegree),

     if verbose,
      disp(['XAB: User-supplied degree less than ' ...
           'minimum expected degree (' int2str(mindegree) ').']);
     end;
     solution = 0; % no solution of given degree
 
    else

     % Upper bound, cf. [1].

     cdeg = max([cdA; cdB], [], 1);
     cdeg(cdeg < 0) = 0; cdeg = [-sort(-cdeg) zeros(1, rA+rB-cA)];
     rdeg = [deg(A, 'row'); deg(B, 'row')];
     rdeg(rdeg < 0) = 0; rdeg = [-sort(-rdeg); zeros(cA-rA-rB, 1)];
     maxdegree = min(sum(cdeg(1:(rA+rB-1))), sum(rdeg(1:(rA+rB-1))));

     if ~isempty(degree) & (degree > maxdegree),
      if verbose,
       disp(['XAB: User-supplied degree greater than ' ...
             'maximum expected degree (' int2str(maxdegree) ').']);
      end;
      maxdegree = degree;
     end;  

    end;

   else % Constant solution.

    mindegree = 0;
    maxdegree = 0;

   end;

   if solution,

    if verbose,
     if ~isempty(degree),
      if degree >= 0,
       disp(['XAB: Seek solution of degree ' int2str(degree) '.']);
      else
       disp(['XAB: Seek solution of degree ' int2str(maxdegree) '.']);
      end;
     elseif maxdegree > 0,
      disp('XAB: Seek minimum degree solution.');
     else
      disp('XAB: Seek constant solution.');
     end;
    end;

    % Binary search for minimum degree
    % ................................
   
    if isempty(degree), degX = mindegree;
    elseif degree < 0,  degX = maxdegree;
    else                degX = degree; end;

    lowerbound = degX;
    upperbound = maxdegree;

    solution = 0; lost = 0; givendegree = ~isempty(degree); stay = 1;

    iteration = 0;

    while stay,

     if iteration,
      degX = ceil((upperbound + lowerbound) / 2);
     end;

     iteration = iteration + 1;

     if verbose,
      fprintf('XAB: Attempt #%3d  Degree %3d (Min %3d / Max %3d) ', ...
       iteration, degX, lowerbound, upperbound);
     end;

     % Coefficients of polynomial matrix B into
     % a constant matrix B of consistent dimensions.
     delta = degA + degX - degB; Bcoef = B.c;
     if delta, RB = [Bcoef(:, :) zeros(rB, cB*delta)];
     else, RB = Bcoef(:, :); end;

     RA = sylv(A, degX);

     % Rank comparison in order to detect solution.

     R = svd(RA);
     rkRA = sum(R > max(size(RA)') * max(R) * tolrank);

     RAB = [RA; RB]; R = svd(RAB);
     rkRAB = sum(R > max(size(RAB)') * max(R) * tolrank);

     if rkRA == rkRAB,

      upperbound = degX; % a solution is found
      oldRA = RA; oldRB = RB; oldrkRA = rkRA; olddegX = degX; % keep it
      solution = 1; lost = 0;

      if verbose,
       fprintf('Yes.\n');
      end;

     else

      lowerbound = degX; lost = 1; % no solution

      if verbose,
       fprintf('No.\n');
      end;

     end;

     if ~solution & (degX == maxdegree - 1), stay = 1; % last chance
     else stay = (upperbound - lowerbound > 1); end;
     if givendegree, givendegree = 0; stay = 0; end; % given degree: only once

    end; % while

    if solution,

     if lost,
      % when necessary, restore minimum degree solution
      RA = oldRA; RB = oldRB; rkRA = oldrkRA; degX = olddegX;
     end;

     % Solve linear system.
     % ....................
 
     [rRA, cRA] = size(RA);

     if rkRA == rRA, % full row rank

      [Q,R] = qr(RA');
      X = pol(RB*Q(:,1:rRA)/R(1:rRA,:)', degX); % unique solution
      X.v = typeA; X.h = Xh;
      
      if verbose,
       disp('XAB: Solve equation with QR factorization.');
      end;

     elseif rkRA == cRA, % full column rank
 
      [Q,R] = qr(RA);
      X = pol((RB/R(1:cRA,:))*Q(:,1:cRA)', degX); % one possible solution
      X.v = typeA; X.h = Xh;
      
      if verbose,
       disp('XAB: Solve equation with QR factorization.');
      end;

     else % rank deficient

      [U,R,V] = svd(RA);
      R = diag(R); R = diag(1./R(1:rkRA));
      X = pol(RB*V(:,1:rkRA)*R*U(:,1:rkRA)', degX); % one possible sol.
      X.v = typeA; X.h = Xh;
      
      if verbose,
       disp('XAB: Solve equation with SVD.');
      end;

     end; % rank
    
     X = pzer(X, tolzero); X.h = Xh;

    end; % solution
   end; % solution
  end; % solution
 end; % method

 % no solution
 if ~solution,
  X = pol(ones(rB, rA)*NaN); X.h = Xh;
 end;

end; % special cases

if solution,
 pprop(X, typeA); X.h = Xh;
end;

% ***********************************
% Null-space computation if required.
% ***********************************

if (nargout == 2) | (solution & ~isempty(columns)),

 if verbose,
  disp('XAB: Null-space computation.');
 end;

 K = null(A.', tol).'; K.h = Xh;

end;

% ************************
% Minimize column degrees.
% ************************

if solution & ~isempty(columns),

 if ~isempty(K),

  if verbose,
   disp('XAB: Degree squeezing.');
  end;

  rkA = size(K, 1); % rank of A

  % Extraction of a non-singular submatrix in K, starting from the
  % right-hand side.  Give priority to the columns specified by the user.

  if size(columns, 2) == 1, columns = columns'; end
  columns = ~~columns;
  mycolumns = columns;
  morecolumns = ~columns;
  mycolumns = mycolumns .* (1:rA);
  mycolumns = fliplr(mycolumns(mycolumns > 0));
  morecolumns = morecolumns .* (1:rA);
  morecolumns = fliplr(morecolumns(morecolumns > 0));
  notfull = 0;

  if verbose,
   disp('XAB: Extract non-singular submatrix in null-space.');
  end;

  rkK = 0; subK = pol(zeros(rkA)); i = 1;
  colind = zeros(1, rkA);
  while (rkK < rkA),
   if ~notfull & (i > length(mycolumns)),
    if verbose,
     disp('XAB: Extra degrees will be minimized.');
    end;
    notfull = 1; mycolumns = morecolumns; i = 1;
   elseif i > length(mycolumns),
    error('XAB: Incorrect rank of the polynomial null-space.');
   end;    
   subK(:,1+rkK) = K(:, mycolumns(i));
   colind(rkA-rkK) = mycolumns(i);
   rkK = rank(subK);
   i = i + 1;
  end;

  % Perform degree reduction by right division on the columns
  % corresponding to the non-singular submatrix in K.

  colind = sort(colind);

  if verbose,
   if sum(columns) > length(colind),
    disp('XAB: Some degrees cannot be minimized.');
   end;
   fprintf('XAB: Indices of minimized degrees = ');
   for i = 1:length(colind),
    fprintf('%d ', colind(i));
   end;
   fprintf('\n');
  end;   

  if verbose,
   disp('XAB: Perform degree reduction through polynomial division.');
  end;

  % null-space for degree squeezing
  D = K(:, colind);

  % division for degree squeezing
  [Q, R] = rdiv(X(:, colind), D, tol);

  if isempty(Q),
   warning(['Polynomial division has no solution. ' ...
            'Degrees cannot be reduced.']);
  else

   X = X - Q*K;
   X = pzer(X, tolzero);

  end;

  if verbose,
   fprintf('XAB: Resulting degrees = ');
   coldeg = deg(X, 'col');
   for i = 1:length(coldeg),
    fprintf('%d ', coldeg(i));
   end;
   fprintf('\n');
  end;   

 elseif verbose,

  disp('Cannot perform degree minimization: no polynomial null-space.');

 end;

end; % minimize

% ****************
% Return solution.
% ****************

if solution,

 if verbose,
  degX = max(0, deg(X));
  disp(['XAB: A solution of degree ' int2str(degX) ' was found.']);
 end;

else

 % No solution.

 if verbose,
  disp('XAB: No polynomial solution was found.');
 end;

end;

%end .. xab

