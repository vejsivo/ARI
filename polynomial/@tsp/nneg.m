function NN = nneg(T)
%NNEG  Nonnegative part of two-sided polynomial
%
% The command
%    NN = NNEG(T)
% for the two-sided polynomial T, returns the polynomial
% with nonnegative powers only.
%
% See also TSP/NPOS, TSP/NEG, TSP/POS.

%        Author: J. Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.
%        $Revision: 3.0 $  $Date: 29-Sep-1999  13:00  $
%                          $Date: 22-May-2000  11:15  $  

Tp = T.p; Tpc = Tp.c;
To = T.o; Ts = T.s; Tr = T.r;

if To<=0,
   NNc = Tpc(:,:,-To+1:end);
else
   NNc = zeros(Ts(1),Ts(2),To+Tr+1);
   NNc(:,:,To+1:end) = Tpc;
end;
NNr = Tr+To;
if NNr>=0,
   NN = pol(NNc(:,:),NNr,'z');
   NN.h = T.h;
else
   NN = pol(zeros(Ts(1),Ts(2)));
end;

%end .. @tsp/nneg
