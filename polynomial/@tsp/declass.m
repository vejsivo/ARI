function [A,C] = declass(T)
%DECLASS  Declass two-sided polynomial (convert to possible lower class)
%            A = DECLASS(T)
%
% The command
%     A = DECLASS(T)
% declasses a two-sided polynomial matrix, i.e.converts its class to
% double (standard Matlab matrix) if possible (if it is empty or
% constant tsp matrix), or converts its class to polynomial in z or
% in z^-1 if possible.
%
% The command
%     [A,C] = DECLASS(T)
% returns also the resulting class in C.
%
% See also TSP/DOUBLE, TSP/POL, TSP/NNEG, TSP/NPOS.

%        Author:  J. Jezek  18-10-99
%        Copyright (c) 1999 by Polyx, Ltd.

if isempty(T) | T.r<0,
   A = zeros(T.s(1),T.s(2));
elseif T.d==0 & T.t==0,
   A = T.p.c(:,:,1);
elseif T.t>=0,
   A = nneg(T);
elseif T.d<=0,
   A = npos(T);
else
   A = T;
end;
C = class(A);

%end .. @tsp/declass
