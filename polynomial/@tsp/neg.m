function N = neg(T)
%NEG   Negative part of two-sided polynomial
%
% The command
%    N = NEG(T)
% for the two-sided polynomial T, returns the polynomial
% with negative powers only.
%
% See also TSP/POS, TSP/NNEG, TSP/NPOS.

%        Author:  J. Jezek  11-8-99
%        Copyright 1999 by Polyx, Ltd.

U = pos(mirror(T));
N = mirror(U);

%end .. @tsp/neg
