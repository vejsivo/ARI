function test = isinf(X)
%ISINF  Test if fraction is infinite
%
% ISINF(X) for fractions X always returns 0.
% This macro exists only for completeness.
%
% See also: ISFINITE, ISNAN.

%      Author: J.Jezek, 22-Sep-2002
%      Copyright(c) 2002 by Polyx, Ltd.

[s1,s2] = size(X);
test = repmat(logical(0),s1,s2);

%end .. @frac/isinf