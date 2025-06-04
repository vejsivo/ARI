function flag = isnan(T)
%ISNAN   Test is two-sided polynomial is Not-a-Number
%
% ISNAN(X) returns an array that contains ones where the 
% entries of the tsp matrix X are NaNs and zeros 
% where they are not.
%
% See also: TSP/ISINF, TSP/ISFINITE.

%   Author: J. Jezek,  11-8-99
%   Copyright 1998 by Polyx, Ltd.

flag = isnan(T.p);

%end .. @tsp/isnan
