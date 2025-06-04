function F = rdf(R,tol)
%RDF    Convert matrix-den fraction to right-den fraction
%
% The command  F = RDF(R)  converts matrix-denominator fraction R
% to right-denominator fraction  F = N*D^-1 . An optional input
% argument TOL may specify a zeroing toleramce to be used instead
% of the standard one.
%
% See also  RDF/RDF, RDF/MDF.

%      Author:  J. Jezek  06-Jan-2000
%      Copyright(c) by Polyx, Ltd.
%      $ Revision $  $Date 21-Apr-2000 $
%                    $Date 29-May-2000 $
%                    $Date 30-Sep-2002 $
%                    $Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Rs1 = R.frac.s(1); Rs2 = R.frac.s(2);
Fd = pol(zeros(Rs2,Rs2));
Fn = pol(zeros(Rs1,Rs2));
for j = 1:Rs2,
   [Fd(j,j),Delta] = plcm(R.frac.den(:,j),[],tol);
   Fn(:,j) = times(R.frac.num(:,j),Delta,tol);
end;
F = rdf(Fn,Fd);

props(F,R.frac.p,R.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   F = coprime(F,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   F = reduce(F,tol);
else
   F = smreduce(F);
end;

%end .. @mdf/rdf
