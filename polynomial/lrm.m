function L = lrm(varargin)
%LRM   Least right multiple
%
% The command
%    L = LRM(N1, N2, .., Nk) 
% computes a least right common multiple L of several polynomial 
% matrices N1, N2, .., Nk (k > 0) that all have the same numbers
% of rows.
%
% Matrices Mi such that L = Ni*Mi may be recovered with Mi = AXB(Ni,L).
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: LLM.

%     Author: D. Henrion,  May 26, 1998.
%     Updated to 3.0 by D. Henrion, August 30, 2000.
%     Copyright 1998-2000 by Polyx, Ltd.
%     Modified by J.Jezek, 24-Jul-2001, case of zero number of matrices

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

% verbose level
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% Handle the case when LRM is called through LLM with
% transposed arguments.

if nargin>0 & isa(varargin{1}, 'cell'),
 varargin = varargin{1};
end;
%nargin = length(varargin);

% Parse input arguments.
% Arrange input matrices into a cell array N = {N1, N2, ..}.

N = {}; rN = [];
tol = [];
degree = -Inf;
n = 0; % number of input matrices
m = 0; % row dimension

if nargin>0,
 for i = 1:nargin,
  arg = varargin{i}; 
  ismat = 0;
  if isa(arg, 'double') % & (i > 1),
   if any([1 1]-size(arg)),
    ismat = 1;
   else
    tol = arg;
   end;
  elseif isa(arg, 'pol'),
   ismat = 1;
  else
   error(['Invalid ',nth(i),' argument.']);
  end;
  if ismat,
   arg = pol(arg);
   if any(any(isnan(arg))) | any(any(isinf(arg))),
    error('Polynomial is not finite.');
   elseif ~isempty(arg),
    n = n + 1;
    if n == 1,
     [m,p] = size(arg);
    else
     [q,p] = size(arg);
     if q ~= m, 
      error('Matrix has not the same number of rows as the other ones.');
     end;
    end;
    cN(n) = p;
    degree = max(degree, deg(arg));
    if rem(n, 2),
     N{n} = arg; 
    else
     N{n} = -arg;
    end;
   end;
  end;
 end;
end;
if n==0,
 error('Number of matrices is zero.');
end;

% Build matrix P = [N1 -N2 0 .. 0; 0 -N2 N3 0 .. 0; 0 0 0 N3 -N4 0 ..];

scN = cumsum(cN);
P = pol(zeros(m*(n-1), scN(n)));
for i = 1:n-1,
 P(1+i*m:(i+1)*m, scN(i)-cN(i)+1:scN(i)) = N{i};
 P(1+i*m:(i+1)*m, scN(i+1)-cN(i+1)+1:scN(i+1)) = N{i+1};
end;

% Zeroing tolerance

Pc = P.coef; me = norm(Pc(:,:), 'inf');
if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else
 if ~isreal(tol) | tol<0 | tol>1,
  error('Invalid tolerance.');
 end;   
 tolzero = tol * me;
end;

% Compute null space of P.

Q = null(P, tol);

% Extract common multiple.

L = pzer(N{1}*Q(1:scN(1), :), tolzero);

%end .. lrm



