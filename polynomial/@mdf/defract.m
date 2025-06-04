function [Q,C] = defract(F)
%DEFRACT     Defract matrix-den fraction 
%             (remove denominators if possible)
%
% The command  Q = DEFRACT(F)  defracts matrix-den fraction F,
% i.e. removes the denominators if they all are eqaul to 1.
%
% The command  [Q,C] = DEFRACT(F)
% returns also the resulting class in C.
%
% See also MDF/COPRIME, MDF/REDUCE, FRAC/POL.

%      Author:  J. Jezek  25-Apr-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 14-Oct-2002 $

global PGLOBAL;

if strcmp(PGLOBAL.COPRIME,'cop'),
   F = coprime(F);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   F = reduce(F);
end;

if all(all(F.frac.den == 1)),
   Q = F.frac.num;
else
   Q = F;
end;
C = class(Q);

%end .. @mdf/defract

