function dtA = det(A,arg1,arg2)
%DET    Determinant of scalar-den fraction
%
% The command
%    DET(A) 
% computes the determinant of a square scalar-den fraction.
% Optional input arguments for method or for tolerance may
% be given like for polynomial matrices.
%
% See also FRAC/RANK.

%      Author:  J.Jezek, 03-Jan-2000
%      Copyright (c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 26-Apr-2000 $
%                    $ Date 06-Nov-2000 $
%                    $ Date 30-Sep-2002 $
%                    $ Date 14-Oct-2002 $

global PGLOBAL;

ni = nargin; tol = []; met = '';
if ni>=2,
   if isa(arg1,'double'), tol = arg1;
   elseif isa(arg1,'char'), met = arg1;
   else error('Invalid 2nd argument.');
   end;
end;
if ni==3,
   if isa(arg2,'double'), tol = arg2;
   elseif isa(arg2,'char'), met = arg2;
   else error('Invalid 3rd argument.');
   end;
end;

num = 0;
eval('num = det(A.frac.num,tol,met);','error(peel(lasterr));');
den = power(A.frac.den,A.frac.s(1),tol);
dtA = sdf(num,den);

if strcmp(A.frac.p,'prop'), props(dtA,'prop',A.frac.tp);
end;
if strcmp(PGLOBAL.COPRIME,'cop'), dtA = coprime(dtA);
end;
if strcmp(PGLOBAL.REDUCE,'red'), dtA = reduce(dtA);
else dtA = smreduce(dtA);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), dtA = defract(dtA);
end;

%end .. @sdf/det
