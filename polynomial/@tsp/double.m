function A = double(T)
%DOUBLE  Convert two-sided polynomial to double (standard Matlab matrix)
%           A = DOUBLE(T)
%
% The command
%     A = DOUBLE(T)
% converts a two-sided polynomial matrix T to double (standard
% Matlab matrix) if possible (if it is empty or constant two-sided
% polynomial matrix). If not possible, error.
%
% See also TSP/DECLASS.

%        Author: J. Jezek   11-10-99
%        Copyright (c) 1999 by Polyx, Ltd.

if isempty(T) | T.r<0,
   A = zeros(T.s(1),T.s(2));
elseif T.d==0 & T.t==0,
   A = T.p.c(:,:,1);
else
   error('Argument is not convertible to double.');
end;

%end .. @tsp/double

