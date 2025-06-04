function [Q,C] = declass(F)
%DECLASS     Declass left-den fraction (convert to possible lower class)
%
% The command  Q = DECLASS(F)  declasses left-den fraction F, i.e.
% converts it to two-sided polynomial (only in 'z' or 'z^-1'),
% polynomial or double (standard MATLAB matrix) if possible.
%
% The command  [Q,C] = DECLASS(F)
% returns also the resulting class in C.
%
% See also LDF/TSP, FRAC/POL, FRAC/DOUBLE.

%      Author:  J. Jezek  10-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      Revision  $ 20-Sep-2001 $
%                $ 14-Oct-2002 $     

if ~isempty(F),
   FF = reduce(coprime(F));
   [Dd,Dlc] = deg(FF.frac.den,'row');
   if strcmp(FF.frac.v,'z') | strcmp(FF.frac.v,'z^-1'),
      if all(Dd==tdeg(FF.frac.den,'row')),
         Q = FF.frac.num;
         if any(Dd>0),
            Dd = repmat(Dd,1,FF.frac.s(2));
            Q = tsp(shift(Q,-Dd,FF.frac.v));
         end;
         Q = Dlc\Q;
         Q = declass(Q);
      else
         Q = F;
      end;
   else
      if all(Dd==0),
         Q = Dlc \ FF.frac.num;
         Q = declass(Q);
      else
         Q = F;
      end;
   end;
else
   Q = zeros(F.frac.s);
end;
C = class(Q);

%end .. @ldf/declass

