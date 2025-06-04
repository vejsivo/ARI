function P = pos(T)
%POS    Positive part of two-sided polynomial
%
% The command
%     P = POS(T)
% for the tsp matrix T, returns the polynomial with
% positive powers only.
%
% See also TSP/NEG, TSP/NPOS, TSP/NNEG.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.
%        $Revision:3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 22-May-2000  11:15:00  $

Tp = T.p; Tpc = Tp.c;
To = T.o; Ts = T.s; Tr = T.r;
Pr = To+Tr;

if Pr<0,
   P = pol(zeros(Ts(1),Ts(2))); return;
end;
Pc = zeros(Ts(1),Ts(2),Pr+1);
if To<=0,
   Pc(:,:,2:end) = Tpc(:,:,-To+2:end);
else
   Pc(:,:,To+1:end) = Tpc;
end;
P = pol(Pc(:,:),Pr,'z');
P.h = T.h;

%end .. @tsp/pos
