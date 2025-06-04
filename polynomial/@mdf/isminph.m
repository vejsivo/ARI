function s = isminph(A,tol)
%ISMINPH  Test if matrix-den fraction is minimum phase
%
% ISMINPH(A) returns 1 if the square matrix-denominator fraction
% is inversely stable and 0 otherwise. See POL/ISSTABLE
% for more information on stability of a polynomial matrix.
%
% See also MDF/ISSTABLE, FRAC/ISMINPH.
%
%      Author: D. Henrion, August 4, 2000.
%      Modified by J. Jezek, July 8, 2001,
%      Copyright 2000 by Polyx, Ltd.

if nargin == 1,
   A = rdf(A);
   eval('s = isminph(A);', 'error(peel(lasterr));');
else
   if ~isempty(tol) & ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
   A = rdf(A,tol);
   eval('s = isminph(A,tol);', 'error(peel(lasterr));');
end;

%end .. @mdf/isminph
