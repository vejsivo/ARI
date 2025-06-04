function D = double(R)
%DOUBLE    Convert fraction to double (standard Matlab matrix)
%
% The command  D = DOUBLE(R)  converts fraction R
% to double (standard Matlab matrix). If not possible, error.
%
% See also FRAC/DECLASS.

%      Author:  J. Jezek  28-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.

eval('D = double(pol(R));', ...
   'error(''Argument is not convertible to double.'');');

%end .. @frac/double
