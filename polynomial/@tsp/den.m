function D = den(F)
%DEN    Denominator of two-sided polynomial
%
% For two-sided polynomial F, the command   D = DEN(F)
% delivers polynomial D, which is denominator of F,
% considered as a fraction in variable 'z'.
%
% This macro exists only for completeness.
%
% See also TSP/NUM, SDF/DEN.

%       Author: J.Jezek, 02-Feb-2002
%       Copyright(c) 2002 by  Polyx, Ltd.
%       $ Revision $  $ Date: 15-Jul-2002 $

if F.o>=0, D = pol(1);
else D = z^-F.o;
end;

%end .. @tsp/den

