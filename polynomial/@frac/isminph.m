function t = isminph(A,tol)
%ISMINPH  Test if fraction is minimum phase
%
% ISMINPH(A) returns 1 if the numerator polynomial matrix of the square
% matrix fraction A (rdf, ldf or sdf) is stable and 0 otherwise.
% See POL/ISSTABLE for more information on stability
% of a polynomial matrix.
%
% See also: FRAC/ISSTABLE, MDF/ISMINPH.

%    Authors: D. Henrion, August 4, 2000.
%    Copyright 2000 by Polyx, Ltd.
%    Modified by J. Jezek, July 6, 2001.
%                          July 25,2002.
%                          Oct  14,2002.

global PGLOBAL;

if nargin == 1,
   tol = PGLOBAL.ZEROING;
else
   if ~isempty(tol) & ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

Acl = class(A);
if strcmp(Acl,'frac'),
   error('Invalid 1st argument.');
end;

t = 0;
eval('t = isstable(A.num,tol);', 'error(peel(lasterr));');

%end .. @frac/isminph
