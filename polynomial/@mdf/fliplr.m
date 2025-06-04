function Q = fliplr(R)
%FLIPLR    Flip matrix-den fraction in left-right direction
%
% The command  Q = FLIPLR(R)  returns R with its rows preserved and
% its columns flipped in the left-right direction.
%
% See also MDF/FLIPUD, MDF/ROT90.

%      Author:  J. Jezek  03-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $     $ 26-Apr-2000 $
%                       $ 30-Sep-2002 $
%                       $ 14-Oct-2002 $

global PGLOBAL;

Q = mdf(fliplr(R.frac.num),fliplr(R.frac.den));

props(Q,R.frac.c,R.frac.tc,R.frac.r,R.frac.p,R.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   Q = coprime(Q);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   Q = reduce(Q);
else
   Q = smreduce(Q);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   Q = defract(Q);
end;

%end .. @mdf/fliplr
