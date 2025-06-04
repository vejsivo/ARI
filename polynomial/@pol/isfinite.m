function test = isfinite(X)
%ISFINITE  Test if polynomial is finite
%
% ISFINITE(X) returns an array that contains ones where the 
% entries of the polynomial matrix X are finite and zeros where 
% they are not.
%
% See also: ISINF, ISNAN.

%   Author: D. Henrion,  May 26, 1998.
%   Copyright 1998 by Polyx, Ltd.
%   $ Revision $  $ Date 19-Jun-2001  P.Zezula  $

d = X.d;
s = size(X);
test = logical(ones(s));
if min(s) > 0, % non-empty
   if d<0, return
   end;   
   for i = 0:d,
      test = test & isfinite(X.c(:,:,i+1));
   end; 
end;

%end .. @pol/isfinite

