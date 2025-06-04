function F = ldf(R,tol)
%LDF    Convert matrix-den fraction to left-den fraction
%
% The command  F = LDF(R)  converts matrix-denominator
% fraction R to left-denominator fraction  F = D^-1*N .
% An optional input argument TOL may specify a zeroing
% tolerance to be used instead of the standard one.
%
% See also LDF/LDF, LDF/MDF.

%      Author:  J. Jezek  06-Jan-2000
%      Copyright(c) by Polyx, Ltd.
%      $ Revision $  $ Date 21-Apr-2000 $
%                    $ Date 29-May-2000 $
%                    $ Date 30-Sep-2002 $
%                    $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Rs1 = R.frac.s(1); Rs2 = R.frac.s(2);
Fd = pol(zeros(Rs1,Rs1));
Fn = pol(zeros(Rs1,Rs2));
for i = 1:Rs1,
   [Fd(i,i),Delta] = plcm(R.frac.den(i,:),[],tol);
   Fn(i,:) = times(R.frac.num(i,:),Delta,tol);
end;
F = ldf(Fd,Fn);

props(F,R.frac.p,R.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   F = coprime(F,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   F = reduce(F,tol);
else
   F = smreduce(F);
end;

%end .. @mdf/ldf
