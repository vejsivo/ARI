function flag = isempty(R)
%ISEMPTY  Test if fraction is empty
%
% ISEMPTY(R) returns 1 if R is an empty fraction
% and 0 otherwise.
%
% See also FRAC/SIZE.

%   Author: J. Jezek, 26-Jan-2000
%   Copyright(c) 2000 by Polyx, Ltd.
%   $ Revision $  $ Date 14-Oct-2002 $

flag = isempty(R.num);

%end .. @frac/isempty
