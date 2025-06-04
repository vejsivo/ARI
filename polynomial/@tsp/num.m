function N = num(F)
%NUM    Numerator of two-sided polynomial
%
% For two-sided polynomial F, the command   N = NUM(F)
% delivers polynomial N, which is numerator of F,
% considered as a fraction in variable 'z'.
%
% This macro exists only for completeness.
%
% See also TSP/DEN, SDF/NUM.

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date: 15-Jul-2002 $

if F.o>=0, N = pol(F);
else N = F.p;
end;

%end .. @tsp/num
