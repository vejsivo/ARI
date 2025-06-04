function flag = isempty(A)
%ISEMPTY  Test if polynomial is empty
%
% ISEMPTY(A) returns 1 if A is an empty polynomial matrix
% and 0 otherwise.
%
% See also  POL/SIZE.

%   Author: D. Henrion,  May 10, 1998.
%   Copyright 1998 by Polyx, Ltd.

flag = (min(size(A))==0);

%end .. @pol/isempty


