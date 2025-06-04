function s = isminph(A,tol)
%ISMINPH   Test if a constant matrix is minimum phase
%
% ISMINPH(A) returns 1 if the square constant matrix A
% is nonsingular.
%
% This macro exists only for completeness.
% See also FRAC/ISMINPH, MDF/ISMINPH.

%      Author: J. Jezek, 08-Jul-2001
%      Copyright(c) 2001 by Polyx, Ltd.
%      $ Revision $  $ Date 11-Sep-2002 $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
elseif ni==1,
   tol = PGLOBAL.ZEROING;
else
   if ~isempty(tol) & ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

eval('s = isminph(sdf(A),tol);', 'error(peel(lasterr));');

%end .. isminph
