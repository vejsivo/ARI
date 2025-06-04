function A = lmispcof(B,tol)
%LMISPCOF  Polynomial spectral co-factorization via LMIs
%
% Given a para-Hermitian polynomial matrix B positive definite on
% the stability boundary, the command
%    A = LMISPCOF(B)
% computes a stable polynomial matrix A such that B = A*A'.
%
% The spectral factor is computed via LMI optimization, see macro LMISPF.

% Author: Didier Henrion, February 8, 2000.
% Updated to 3.0 by Didier Henrion, September 1, 2000.
% Copyright 2000 Polyx, Ltd.
% Modified by Jan Jezek, August 2001, arg checking

if nargin == 0,
   error('Not enough input arguments.');
elseif nargin == 1,
   eval('A = lmispf(B,[],''cofactor'');', ...
      'error(peel(lasterr));');
elseif nargin == 2,
   eval('A = lmispf(B,tol,''cofactor'');', ...
      'error(peel(lasterr));');
end;

%end .. lmispcof
