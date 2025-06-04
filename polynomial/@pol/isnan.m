function test = isnan(X)
%ISNAN  Test if polynomial is Not-a-Number
%
% ISNAN(X) returns an array that contains ones where the 
% entries of the polynomial matrix X are NaNs and zeros 
% where they are not.
%
% See also: POL/ISINF, POL/ISFINITE.

%   Author: D. Henrion,  May 26, 1998.
%   Copyright 1998 by Polyx, Ltd.
%   $ Revision 3.0 $  $ Date 20-Jun-2001  J.Jezek  $

d = X.d;
s = size(X);
test = logical(zeros(s));
if min(s) > 0, % non-empty
   if d<0, return
   end;
   for i = 0:d,
      test = test | isnan(X.c(:,:,i+1));
   end; 
end;

%end .. @pol/isnan

