function t = isstable(A,tol)
%ISSTABLE  Test if matrix-den fraction is stable
%
% ISSTABLE(A) returns 1 if the denominator polynomials in all entries
% of matrix of A are stable and 0 otherwise. See POL/ISHURWITZ,
% POL/ISSCHUR for more information on stability of a polynomial.
%
% See also: FRAC/ISSTABLE, FRAC/ISMINPH.

%    Authors: D. Henrion, October 9, 2000.
%    Modified by J. Jezek, July 6, 2001.
%    Copyright 2000 by Polyx, Ltd.

global PGLOBAL;

if nargin == 1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

D = A.frac.d;
[sizeD1,sizeD2] = size(D);
t = logical(1);
if sizeD1==0 | sizeD2==0,
   return;
end;

for i = 1:sizeD1,
   for j = 1:sizeD2,
      sij = isstable(D(i,j),tol);
      if ~sij,
         t = logical(0); return;
      end;
   end;
end;

%end .. @mdf/isstable








