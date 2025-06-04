function test = isnan(X)
%ISNAN  Test if fraction is Not-a-Number
%
% ISNAN(X) for fractions X always returns 0.
% This macro exists only for completeness.
%
% See also: ISFINITE, ISINF.

%      Author: J.Jezek, 22-Sep-2002
%      Copyright(c) 2002 by Polyx, Ltd.

[s1,s2] = size(X);
test = repmat(logical(0),s1,s2);

%end .. @frac/isnan