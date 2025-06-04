function M = sdf(F,tol)
%SDF   Convert left=den fraction to scalar-den fraction
%
% The command  M = SDF(F)  converts left-denominator fraction
% F to scalar-denominator fraction M. The denominator of M is
% the determinant of the denominator of F.
%
% An optional input argumenr TOL may specify a zeroing tolerance
% to be used instead the standard one.
%
% See also SDF/SDF, SDF/LDF.

%       Author:  J. Jezek  24-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 02-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

AdjFd = 0;
eval('[AdjFd,Det] = adj(F.frac.den,tol);','error(peel(lasterr));');
N = mtimes(AdjFd,F.frac.num,tol);
M = sdf(N,Det);

props(M,F.frac.p,F.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   M = coprime(M,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   M = reduce(M);
else
   M = smreduce(M);
end;

%end .. @ldf/sdf


