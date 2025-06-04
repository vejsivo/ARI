function flag = isempty(T)
%ISEMPTY  Test if two-sided polynomial is empty
%
% ISEMPTY(A) returns 1 if A is an empty tsp matrix
% and 0 otherwise.
%
% See also  TSP/SIZE.

%   Author: J. Jezek, 11-8-99
%   Copyright 1998 by Polyx, Ltd.

flag = isempty(T.p);

%end .. @tsp/isempty
