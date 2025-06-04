function Q = reshape(P,M,N)
%RESHAPE        Change size of matrix-den fraction
%
% For matrix-den fraction P, the command Q = RESHAPE(P,M,N)  returns
% a M-by-N matrix-den fraction whose elements are taken columnwise
% from P. An error results if P has not M*N elements.
%
% Q = RESHAPE(P,[M N])  is the same thing.

%      Author:  J. Jezek  02-Feb-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 26-Apr-2000 $
%                    $ Date 06-Nov-2000 $
%                    $ Date 30-Sep-2002 $
%                    $ Date 14-Oct-2002 $

global PGLOBAL;

eval('P = mdf(P);','error(peel(lasterr));');

ni = nargin;
if ni==1,
   error('Not enough input arguments.');
end;
if ~isa(M,'double'),
   error('Invalid 2nd argument.');
end;

Qn = 0;
if ni==2,
   eval('Qn = reshape(P.frac.num,M);','error(peel(lasterr));');
   Qd = reshape(P.frac.den,M);
else
   if ~isa(N,'double'),
      error('Invalid 3rd argument.');
   end;
   eval('Qn = reshape(P.frac.num,M,N);','error(peel(lasterr));')
   Qd = reshape(P.frac.den,M,N);
end;
Q = mdf(Qn,Qd);

props(Q,P.frac.p,P.frac.tp,P.frac.r,P.frac.c,P.frac.tc);
if strcmp(PGLOBAL.COPRIME,'prop'), Q = coprime(Q);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Q = reduce(Q);
else Q = smreduce(Q);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Q = defract(Q);
end;

%end .. @mdf/reshape
