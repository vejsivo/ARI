function R = mdf(F,tol)
%MDF    Convert right-den fraction to matrix-den fraction
%
% For right-denominator fraction F, the command  R = MDF(F)
% returns matrix-denominator fraction R, equal to F.
%
% An optional input argument TOL may specify a zeroing tolerance to be 
% used istead the standard one.
%
% See also MDF/MDF, MDF/RDF.

%        Author:  J. Jezek  06-Jan-2000
%        Copyright(c) by Polyx, Ltd.
%        $ Revision $  $ Date 21-Apr-2000 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

[Dadj,Ddet] = adj(F.frac.den,tol);
R = mdf(mtimes(F.frac.num,Dadj,tol),Ddet);

props(R,F.frac.p,F.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   R = coprime(R,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   R = reduce(R);
else
   R = smreduce(R);
end;

%end .. @rdf/mdf
