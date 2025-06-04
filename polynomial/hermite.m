function [T,U,rowind] = hermite(A, varargin)
%HERMITE  Hermite form of a polynomial matrix
%
% The commmands
%    T = HERMITE(A)
%    T = HERMITE(A,'col') 
% compute the column Hermite form of A. The Hermite form is unique. 
% It features monic leading entries and reduced degree non-leading entries.
%    
% The commmand
%    [T,U,ROWIND] = HERMITE(A) 
% additionally provides a unimodular reduction matrix U such that AU = T,
% and a row vector ROWIND such that ROWIND(I) is the row index of the 
% uppermost nonzero entry in the Ith column of T. Note that ROWIND(I) = 0 
% if the I-th column of T is zero.
%
% Similarly, 
%   [T,U,COLIND] = HERMITE(A,'row') 
% computes the row Hermite form of A, the unimodular reduction matrix U such 
% that UA = T, and a column vector COLIND such that COLIND(I) is the column 
% index of the leftmost nonzero entry in the Ith row of T.
%
% The macro HERMITE is partly based on the macro TRI. As a result, the optional
% input arguments for TRI may be specified through HERMITE.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: TRI.

%    Author: D. Henrion, September 21, 1998.
%    Updated to 3.0 by D. Henrion, August 30, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin<1,
   error('Not enough input arguments.');
end;
eval('A = pol(A);', 'error(peel(lasterr));');

if any(any(isnan(A))) | any(any(isinf(A))),
 error('Polynomial is not finite.');
end;

% verbose level
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
    else
     error('Invalid command option.');
    end;     
   elseif isa(arg, 'double'), % matrix or scalar
    if ~any(size(arg) - [1 1]), % scalar
     if floor(arg) ~= arg, % non integer = tolerance
      tol = arg;
     else % integer = degree
      degree = arg;
     end;
    else % matrix argument
     error(['Invalid ',nth(i+1),' argument.']);
    end;
   else
    error(['Invalid ',nth(i+1),' argument.'])    
   end;
  end;
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

% zeroing tolerance
Ac = A.coef;
me = norm(Ac(:,:), 'inf');
if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else
 tolzero = tol * me;
end;

if isempty(A),

 % Empty matrix.
 T = pol(zeros(cA, 0));
 U = pol(zeros(rA, 0));
 rowind = [];
 
elseif dA < 0,

 % Zero matrix.
 T = pol(zeros(rA, cA));
 U = pol(eye(cA));
 rowind = [];

else,

 % ****************************************************************
 % First step: call TRI to obtain a non-canonical staircase form.
 % ****************************************************************

 [T,U,rowind] = tri(A,method,degree,tol);

 % ************************************************************
 % Second step: transform non-canonical form into Hermite form.
 % ************************************************************

 if verbose,
  disp('HERMITE: Transformation into Hermite form.');
 end;

 V = pol(eye(cA)); % reduction matrix

 % Reduce off diagonal entry degrees with RDIV.
 % Perform zeroing.

 col = 1;
 RI = rowind(2:end);
 for row = RI(:)',
  if row > 0,

   if verbose,
    disp(['HERMITE: Reduction of column #' int2str(col) '.']);
   end;

   col = col + 1;
   N = T(row, 1:col-1); % non-leading entries
   D = T(row, col); % leading entry

   [Q, R] = rdiv(N, D, tol);
   if any(any(isnan(Q))) | any(any(isnan(R))),
    error('The reduction failed. Try to modify the tolerance.');
   end;

   % refresh non-leading entries in row
   T(row, 1:col-1) = pzer(R, tolzero);

   if row < rA,
    % refresh remaining rows
    T(row+1:rA, 1:col-1) = pzer(T(row+1:rA, 1:col-1) - ...
     T(row+1:rA, col) * Q, tolzero);
   end;

   % refresh reduction matrix
   V(col, 1:col-1) = -pzer(Q, tolzero);

  end;
 end;

 % Make leading entries monic. 

 col = 0;
 for row = rowind(:)',
  if row > 0,
   col = col + 1;
   coef = lcoef(T(row, col));
    
   if abs(coef) < tolzero,
    warning(['Hermite form features a small leading coefficient. ' ...
             'Results may be inaccurate.']);
   end;

   % refresh Hermite form
   T(row:rA, col) = T(row:rA, col) / coef;

   % refresh reduction matrix
   V(:, col) = V(:, col) / coef;

  end;
 end;

 % Update reduction matrix.

 U = pzer(U*V, tolzero);

 dU = deg(U);
 dT = deg(T);

 if isinf(dU) | isinf(dT),
  disp('The reduction matrix or Hermite form is zero. The reduction failed.');
  error('Try to modify the tolerance.');
 end;

 if verbose,
  disp(['HERMITE: Reduction matrix of degree ' int2str(deg(U)) '.']);
  disp(['HERMITE: Hermite form of degree ' int2str(deg(T)) '.']);
 end;

end;

% ************
% LAST TUNINGS
% ************

if strcmp(form, 'row'),
 U = U.';
 T = T.';
 rowind = rowind';
end;

%end .. hermite

