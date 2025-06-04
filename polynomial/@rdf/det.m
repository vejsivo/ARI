function dtA = det(A,arg1,arg2)
%DET    Determinant of right-den fraction
%
% The command
%    DET(A) 
% computes the determinant of a square right-denominator fraction.
% Optional input arguments for method or for tolerance may be given
% like for polynomial matrices.
%
% See also FRAC/RANK.

%      Author:  J.Jezek, 01-Dec-1999
%      Copyright (c) 1998 by Polyx, Ltd.
%      $ Revision $  $ Date 25-Apr-2000 $
%                    $ Date 01-Nov-2000 $
%                    $ Date 30-Sep-2002 $
%                    $ Date 14-Oct-2002 $

global PGLOBAL;

ni = nargin;
if ni>=2,
   if ~isa(arg1,'double') & ~isa(arg1,'char'),
      error('Invalid 2nd argument.');
   end;
end;
if ni==3,
   if ~isa(arg2,'double') & ~isa(arg2,'char'),
      error('Invalid 3rd argument.');
   end;
end;

num = 0; den = 0;
if ni==1,
   eval('num = det(A.frac.num); den = det(A.frac.den);','error(peel(lasterr));');
elseif ni==2,
   eval('num = det(A.frac.num,arg1); den = det(A.frac.den,arg1);',...
        'error(peel(lasterr));');
else
   eval('num = det(A.frac.num,arg1,arg2); den = det(A.frac.den,arg1,arg2);',...
        'error(peel(lasterr));');
end;

dtA = rdf(num,den);

if strcmp(A.frac.p,'prop'), props(dtA,'prop',A.frac.tp);
end;
if strcmp(PGLOBAL.COPRIME,'cop'), dtA = coprime(dtA);
end;
if strcmp(PGLOBAL.REDUCE,'red'), dtA = reduce(dtA);
else dtA = smreduce(dtA);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), dtA = defract(dtA);
end;

%end .. @rdf/det
