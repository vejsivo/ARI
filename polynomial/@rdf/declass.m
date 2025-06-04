function [Q,C] = declass(F)
%DECLASS     Declass right-den fraction (convert to possible lower class)
%
% The command  Q = DECLASS(F)  declasses right-den fraction F, i.e.
% converts it to two-sided polynomial (only in 'z' or 'z^-1'),
% polynomial or double (standard MATLAB matrix) if possible.
%
% The command  [Q,C] = DECLASS(F)
% returns also the resulting class in C.

% See also RDF/TSP, FRAC/POL, FRAC/DOUBLE.

%      Author:  J. Jezek  10-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      Revision  $ 20-Sep-2001 $
%                $ 14-Oct-2002 $

if ~isempty(F),
   FF = reduce(coprime(F));
   [Dd,Dlc] = deg(FF.frac.den,'col');
   if strcmp(FF.frac.v,'z') | strcmp(FF.frac.v,'z^-1'),
      if all(Dd==tdeg(FF.frac.den,'col')),
         Q = FF.frac.num;
         if any(Dd>0),
            Dd = repmat(Dd,FF.frac.s(1),1);
            Q = tsp(shift(Q,-Dd,FF.frac.v));
         end;
         Q = Q/Dlc;
         Q = declass(Q);
      else
         Q = F;
      end;
   else
      if all(Dd==0),
         Q = FF.frac.num / Dlc;
         Q = declass(Q);
      else
         Q = F;
      end;
   end;
else
   Q = zeros(F.frac.s);
end;
C = class(Q);

%end .. @rdf/declass

