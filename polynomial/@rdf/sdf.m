function M = sdf(F,tol)
%SDF   Convert right-den fraction to scalar-den fraction
%
% The command  M = SDF(F)  converts right-denominator fraction
% F to scalar-denominator fraction M. The denominator of M is
% the determinant of the denominator of F.
%
% An optional input argumenr TOL may specify a zeroing tolerance
% to be used instead the standard one.
%
% See also SDF/SDF, SDF/RDF.

%       Author:  J. Jezek  24-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 21-Apr-2000 $
%                     $ Date 01-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Revision $  $ Date 14-Oct-2002 $


global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

AdjFd = 0; Det = 0;
eval('[AdjFd,Det] = adj(F.frac.den,tol);','error(peel(lasterr));');
N = mtimes(F.frac.num,AdjFd,tol);
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

%end .. @rdf/sdf


