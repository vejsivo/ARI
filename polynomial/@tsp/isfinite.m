function flag = isfinite(T)
%ISFINITE   Test if two-sided polynomial is finite
%
% ISFINITE(X) returns an array that contains ones where the 
% entries of the tsp matrix X are finite and zeros where 
% they are not.
%
% See also: TSP/ISINF, TSP/ISNAN.

%   Author: J. Jezek,  11-8-99
%   Copyright 1998 by Polyx, Ltd.

flag = isfinite(T.p);

%end .. @tsp/isfinite
