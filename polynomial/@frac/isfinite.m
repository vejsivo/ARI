function test = isfinite(X)
%ISFINITE  Test if fraction is finite
%
% ISFINITE(X) for fraction X always returns 1
% of corresponding dimension.
% This macro exists only for completeness.
%
% See also: ISINF, ISNAN.

%      Author: J.Jezek, 22-Sep-2002
%      Copyright(c) 2002 by Polyx, Ltd.

[s1,s2] = size(X);
test = repmat(logical(1),s1,s2);

%end .. @frac/isfinite