function [Q,C] = defract(F)
%DEFRACT     Defract left-den fraction 
%             (remove denominator if possible)
%
% The command  Q = DEFRACT(F)  defracts left-den fraction F,
% i.e. removes the denominator if it is unity matrix.
%
% The command  [Q,C] = DEFRACT(F)
% returns also the resulting class in C.
%
% See also LDF/COPRIME, LDF/REDUCE, FRAC/POL.

%      Author:  J. Jezek  25-Apr-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ 14-Oct-2002 $

global PGLOBAL;

if strcmp(PGLOBAL.COPRIME,'cop'),
   F = coprime(F);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   F = reduce(F);
end;

if all(all(F.frac.den == pol(eye(F.frac.s(1))))),
   Q = F.frac.num;
else
   Q = F;
end;
C = class(Q);

%end .. @ldf/defract

