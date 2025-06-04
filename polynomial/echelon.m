function [E,U,indices] = echelon(D,varargin)
%ECHELON  Echelon (or Popov) form of a polynomial matrix
%
% Given a square and column-reduced polynomial matrix D of 
% dimension N, the commands
%    [E,U,ROWIND] = ECHELON(D)
%    [E,U,ROWIND] = ECHELON(D,'col')
% compute the column polynomial echelon form E of D with
% - its column degrees arranged in ascending order
% - the degree of every pivot element
%   - larger than the degree of every lower element in the same column,
%   - larger than or equal to the degree of every upper element in the
%     same column,
%   - larger than the degree of every other element in the same row.
% Clearly E is both column reduced and row reduced.
%
% The unimodular matrix U satisfies
%    DU = E
% while
%    ROWIND(1:N,1) contains the pivot (row) indices of E,
%    ROWIND(1:N,2) contains the column degrees of E,
%    ROWIND(1:N,3) contains the row degrees of E.
%
% Similarly, given a square and column-reduced polynomial matrix D of 
% dimension N, the command
%    [E,U,COLIND] = ECHELON(D,'row')
% produces E in column polynomial echelon form, satisfying the 
% properties as listed with 'column' replaced by 'row', and vice-versa.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TRI, HERMITE.

%   Authors: D. Henrion, September 21, 1998.
%   Copyright 1998 by Polyx, Ltd.
%   $ Revision 3.0 $  $ Date 30-May-2000  J.Jezek  $
%                     $ Date 28-Jan-2001  J.Jezek  $

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin<1
   error('Not enough input arguments.');
end;
eval('D = pol(D);', 'error(peel(lasterr));');

if any(any(isnan(D))),
 error('Polnomial is not finite.');
elseif any(any(isinf(D))),
 error('Polynomial is not finite.');
end;

% verbose flag
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% Default options.

tol = [];
form = 'col';

% Parse optional input arguments.

invalid = 0;
if nargin > 1,
 for i = 1:length(varargin),
  arg = varargin{i};
  if ~isempty(arg),
   if isa(arg, 'char'), % valid strings
    if strcmp(arg, 'col'),
     form = 'col';
    elseif strcmp(arg, 'row'),
     form = 'row';
    else,
     error('Invalid command option.');
    end;     
   elseif isa(arg, 'double'), % matrix or scalar
    if ~any(size(arg) - [1 1]), % scalar = tolerance
     tol = arg;
    else, % matrix argument
     error('Invalid matrix argument.');
    end;
   else,
    error('Invalid input argument.');
   end;
  end;
 end;
end;

if isempty(tol),
 tol = PGLOBAL.ZEROING;
end;

[n, cD] = size(D);
degD = deg(D);
var = D.var;

if n ~= cD,
  error('Matrix is not square.');
elseif rank(D, tol) < n,
  error('Matrix is singular.');
end;

Dcol = lcoef(D, 'col');
Drow = lcoef(D, 'row');

if (rank(Drow) < n) & strcmp(form, 'row')
  error('Matrix is not row reduced.')
elseif (rank(Dcol) < n) & strcmp(form, 'col')
  error('Matrix is not column reduced.')
end;

if strcmp(form,'col'),
  D = D.';
end;

% Resultants of increasing degree are built until all pivot elements
% are found (i.e. until there is no more "primary" dependent rows).

shft = 0;
notfound = 1:n; % absent pivot indices

while any(notfound),

  degS = 0; k = 1;
  notfound = 1:n; % absent pivot indices
  pdr = zeros(n, 1); % primary dependent rows indices
  notproper = 0; % properness of G(s) = D^-1(s)

  while ~notproper & any(notfound),

    % Sylvester resultant matrix

    R = sylv([D; -eye(n)], degS);

    % Right null-space computation

    [S, irow] = cef(R, tol);
    S = nullref(S, irow);

    % Indices of dependent rows

    rS = size(S, 1); rR = size(R, 1);
    dr = [1:rR]'; dr(irow) = 0; dr = dr(dr > 0);

    degS = degS + 1;
    pivot = rem(dr-1, 2*n);

    if any(pivot < n), % linearly dependent row in the upper part ..

       notproper = 1; % .. so G(s) is not proper
       D = shift(D, 1); % compute G(s) = d^(-shft)D^-1(s)
       shft = shft + 1;

    else

       pivot = pivot-n+1; % pivot indices
       for i = 1:rS,
         j = find(~(notfound - pivot(i)));
         if ~isempty(j), % primary dependent row
           notfound(j(1)) = 0;
           pdr(k) = i; k = k + 1;
         end;
       end;

    end;
  end;
end;

% selection of the primary dependent rows
S = pol(S(pdr, :), degS - 1, var);

% pivot indices, row degrees and column degrees
pivot = pivot(pdr);
rdeg = fix((dr(pdr)-1)/(2*n)) - shft;
[void, ind] = sort(pivot); cdeg = rdeg(ind);
indices = [pivot rdeg cdeg];

U = pzer(S(1:n, 1:n),tol);  % argument TOL added by J.Jezek 28-Jan-2001
E = shift(pzer(S(1:n, n+1:2*n)), -shft);

U.h = D.h; E.h = D.h;

if strcmp(form, 'col'),
 U = U.'; E = E.';
end;

%end .. echelon
