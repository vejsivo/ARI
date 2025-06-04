function A = spcof(B,varargin)
%SPCOF  Two-sided polynomial spectral co-factorization
%
% If B is a para-Hermitian two-sided polynomial matrix,
% that is,   B'(z^-1) = B(z) , positive definite on the unit circle,
% the command
%    A = SPCOF(B)
% returns a stable polynomial A(z), satisfying
%    A(z^-1)*A'(z) = B(z) .
% 
% For more information, see TSP/SPF.

%        Author:  J. Jezek, 11-8-1999
%        Copyright (c) 1999 by Polyx, Ltd.

eval('A = spf(B.'',varargin{:});','error(peel(lasterr));');
A = A.';

%end .. @tsp/spcof


