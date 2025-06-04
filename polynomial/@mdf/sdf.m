function M = sdf(R,tol)
%SDF  Convert matrix-den fraction to scalar-den fraction
%
% The command  M = SDF(R)  converts matrix-denominator fraction
% R to scalar-denominator fraction M. The denominator of M is
% the least common multiple of all denominators of R.
%
% An optional input argument TOL may specify a zeroing tolerance
% to be used instead the standard one.
%
% See also SDF/SDF, SDF/MDF.

%       Author:  J. Jezek  24-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 21-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

DD = 0; Delta = 0;
eval('[DD,Delta] = plcm(R.frac.den,[],tol);','error(peel(lasterr));');
NN = times(R.frac.num,Delta,tol);
[DD,Delta] = plcm(DD,[],tol);
Delta = repmat(Delta,R.frac.s(1),1);
NN = times(NN,Delta,tol);
M = sdf(NN,DD);

props(M,R.frac.p,R.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   M = coprime(M,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   M = reduce(M);
else
   M = smreduce(M);
end;

%end .. @mdf/sdf
