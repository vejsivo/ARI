function [A,C] = declass(P)
%DECLASS  Declass polynomial (convert to possible lower class)
%            A = DECLASS(P)
%
% The command
%     A = DECLASS(P)
% declasses a polynomial matrix P, i.e. changes its class to double
% (standard Matlab matrix) if possible (if it is empty or constant
% polynomial matrix).
%
% The command
%     [A,C] = DECLASS(P)
% returns also the resulting class in C.
%
% See also POL/DOUBLE.

%         Author: J. Jezek  18-10-99
%         Copyright (c) 1999 by Polyx, Ltd.

if isempty(P) | P.d<0,
   A = zeros(P.s(1),P.s(2));
elseif P.d==0,
   A = P.c(:,:,1);
else
   A = P;
end;
C = class(A);

%end .. @pol/declass

