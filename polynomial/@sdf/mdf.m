function R = mdf(M,tol)
%MDF  Convert scalar-den fraction to matrix-den fraction
%
% The command  R = MDF(M)  converts scalar-denominator
% fraction M to matrix-denominator fraction R.
%
% An optional input argument TOL may specify the zeroing
% tolerance to be used instead of the standard one.
%
% See also MDF/MDF, MDF/SDF.

%       Author:  J. Jezek  24-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $Date 21-Apr-2000 $
%                     $Date 30-Sep-2002 $
%                     $Date 14-Oct-2002 $

global PGLOBAL;

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING;
end;

R = mdf(M.frac.num,M.frac.den);

props(R,M.frac.p,M.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   R = coprime(R,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   R = reduce(R,tol);
else
   R = smreduce(R);
end;

%end .. @sdf/mdf
