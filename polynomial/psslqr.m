function L = psslqr(F, G, H, J, tol)
%PSSLQR  Polynomial approach to linear-quadratic regulator design
%         for state-space system
%
% Given a linear system
%   .
%   x = F*x + G*u
%
% where F is an NxN constant matrix and G is an NxM constant matrix, and
% a regulated variable
%
%   z = H*x + J*u
%
% where H is a PxN constant matrix and J is a PxM constant matrix,
% the command
%
%   L = PSSLQR(F, G, H, J)
%
% returns a constant matrix L such that the control function u = -L*x
% minimizes the L2-norm of z for every initial state x(0).
% It is assumed that
%
%   J'*H = 0, J'*J = I.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%     Author: D. Henrion, September 12, 1999.
%     Copyright 1999 by Polyx, Ltd.
%     Modified by J. Jezek, August 2001, arg checking

global PGLOBAL
eval('PGLOBAL.VERBOSE;', 'painit;');

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% ---------------------
% Parse input arguments
% ---------------------

if nargin < 4,
  error('Not enough input arguments.');
end;

if nargin < 5 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isnumeric(tol) | length(tol)~=1 | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

if ~isnumeric(F) | ~isnumeric(G) | ~isnumeric(H) | ~isnumeric(J) | ...
      ndims(F)>2 | ndims(G)>2 | ndims(H)>2 | ndims(J)>2,
   error('Invalid 1st - 4th arguments.');
end;

n = size(F, 1);
m = size(G, 2);
p = size(H, 1);
if (size(F, 2) ~= n) | (size(G, 1) ~= n) | (size(H, 2) ~= n) | ...
   (size(J, 1) ~= p) | (size(J, 2) ~= m),
  error('Matrices of inconsistent dimensions.');
end;

if norm(J'*H) > tol,
 error('Assumption J''*H = 0 does not hold.');
end;

if norm(J'*J-eye(m)) > tol,
 error('Assumption J''*J = I does not hold.');
end;

% ------------------------
% Computation of right MFD
% ------------------------

if verbose,
 disp('PSSLQR: Compute MFD');
end;

[BR, AR] = ss2rmf(F, G, eye(n), zeros(n,m), tol);

% --------------------------------------------------------
% Compute right hand-side matrix by spectral factorization
% --------------------------------------------------------

if verbose,
 disp('PSSLQR: Build right hand-side polynomial matrix by spectral factorization.');
end;

CR = spf(BR'*H'*H*BR+AR'*AR, tol);

% ---------------------------------------------------
% Solve the polynomial matrix equation XL*AR+YL*BR=CR
% for constant matrices XL, YL
% ---------------------------------------------------

if verbose,
 disp('PSSLQR: Solve Diophantine equation.');
end;

[XL, YL] = xaybc(AR, BR, CR, 0, tol);

% ----------------------------
% Build feedback gain matrix L
% ----------------------------

L = XL\YL;

%end .. psslqr


