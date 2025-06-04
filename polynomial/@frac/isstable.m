function t = isstable(A,tol)
%ISSTABLE  Test if fraction is stable
%
% For fraction A (right-den fraction, left-den fraction
% or scalar-den fraction), ISSTABLE(A) returns 1 if the denominator
% polynomial matrix of A is stable and 0 otherwise.
%
% See POL/ISSTABLE for more information
% on stability of a polynomial matrix.
%
% See also: MDF/ISSTABLE, FRAC/ISMINPH.

%    Authors: D. Henrion, August 4, 2000.
%    Last modified by D. Henrion, October 20, 2000.
%                  by J. Jezek, July 06, 2001.
%                               July 25, 2002.
%                               Oct  14, 2002.
%    Copyright 2000 by Polyx, Ltd.

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

% stability test on the denominator
t = 0;
eval('t = isstable(A.den,tol);', 'error(peel(lasterr));');

%end .. @frac/isstable
