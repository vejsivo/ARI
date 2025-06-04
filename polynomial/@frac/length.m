function l = length(P)
%LENGTH  Length of fraction
%   
% L = LENGTH(P) returns the length of fraction P.
% It is equivalent to MAX(SIZE(P)) for a nonempty matrix
% and 0 for an empty matrix.

%       Author:  J. Jezek  08-Feb-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 14-Oct-2002 $

l = length(P.num);

%end .. @frac/length

