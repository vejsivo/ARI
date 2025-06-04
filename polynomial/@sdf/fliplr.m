function Q = fliplr(R)
%FLIPLR    Flip scalar-den fraction in left-right direction
%
% The command  Q = FLIPLR(R)  returns R with its rows preserved and
% its columns flipped in the left-right direction.
%
% See also SDF/FLIPUD, SDF/ROT90.

%      Author:  J. Jezek  26-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 26-Apr-2000 $
%                    $ Date 30-Sep-2002 $
%                    $ Date 14-Oct-2002 $

global PGLOBAL;

Q = sdf(fliplr(R.frac.num),R.frac.den);

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

%end .. @sdf/fliplr
