function flag = isinf(T)
%ISINF   Test if two-sided polynomial is infinite
%
% ISINF(X) returns an array that contains ones where the 
% entries of the tsp matrix X are +Inf or -Inf and
% zeros where they are not.
%
% See also: TSP/ISFINITE, TSP/ISNAN.

%   Author: J. Jezek, 11-8-99
%   Copyright 1998 by Polyx, Ltd.

flag = isinf(T.p);

%end .. @tsp/isinf
