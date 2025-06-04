function NP = npos(T)
%NPOS  Nonpositive part of two-sided polynomial
%
% The command
%    NP = npos(T)
% for the two-sided polynomial T, returns the polynomial
% with nonpositive powers only.
%
% See also TSP/NNEG, TSP/POS, TSP/NEG.

%        Author:  J.Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.

U = nneg(mirror(T));
NP = mirror(U);

%end .. @tsp/npos
