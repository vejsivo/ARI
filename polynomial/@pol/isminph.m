function s = isminph(A,tol)
%ISMINPH   Test whether a polynomial is minimum phase
%
% ISMINPH(A) returns i if the square polynomial matrix A,
% considered as a fraction, is minimum phase.
%
% This macro exists only for completeness.
% See also FRAC/ISMINPH, MDF/ISMINPH.

%      Author: J. Jezek, 08-Jul-2001
%      Copyright(c) 2001 by Polyx, Ltd.

global PGLOBAL;

if nargin == 1,
   tol = PGLOBAL.ZEROING;
else
   if ~isempty(tol) & ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

eval('s = isminph(sdf(A),tol);', 'error(peel(lasterr));');

%end .. @pol/isminph

