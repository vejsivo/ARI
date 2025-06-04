function F = rdf(M,tol)
%RDF  Convert scalar-den fraction to right-den fraction
%
% The command  F = RDF(M)  converts scalar-denominator
% fraction M to right=ddenominator fraction F.
%
% An optional input argument TOL may specify the zeroing
% tolerance to be used instead the standard one.
%
% See also RFD/RDF, RDF/SDF.

%       Author:  J. Jezek  24-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 21-Apr-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING;
end;

F = rdf(M.frac.num,M.frac.den);

props(F,M.frac.p,M.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   F = coprime(F,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   F = reduce(F,tol);
else
   F = smreduce(F);
end;

%end .. @sdf/rdf
