function [Q,C] = declass(R)
%DECLASS     Declass scalar-den fraction (convert to possible lower class)
%
% The command  Q = DECLASS(R)  declasses scalar-den fraction R
% i.e converts it to two-sided polynomial (only in 'z' or 'z^-1'),
% polynomial or double (standard Matlab matrix) if possible.
%
% The command  [Q,C] = DECLASS(R)
% returns also the resulting class in C.

% See also FRAC/TSP, FRAC/POL, FRAC/DOUBLE.

%      Author:  J. Jezek  07-Nov-2002
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 23-Nov-2002  bug $

if ~isempty(R),
   RR = coprime(R);
   [Dd,Dlc] = deg(RR.frac.den);
   if strcmp(RR.frac.v,'z') | strcmp(RR.frac.v,'z^-1'),
      if Dd==tdeg(RR.frac.den),
         Q = RR.frac.num./Dlc;
         if Dd>0,
            Q = tsp(shift(Q,-Dd,RR.frac.v));
         end;
         Q = declass(Q);
      else
         Q = R;
      end;
   else
      if Dd==0,
         Q = RR.frac.num./Dlc;
         Q = declass(Q);
      else
         Q = R;
      end;
   end;
else
   Q = zeros(R.frac.s);
end;
C = class(Q);

%end .. @sdf/declass
