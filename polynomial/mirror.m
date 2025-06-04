function B = mirror(A)
%MIRROR    Mirror image of  constant
%
%  B = MIRROR(A)  for constant A (i.e. the standard
% Matlab matrix, returns  B = A.
%
% This macro exists only for completeness.
% See also POL/MIRROR, FRAC/MIRROR.

%     Author:  J. Jezek, 24-Feb-2003
%     Copyright(c) 2003 by Polyx, Ltd.

if ~isa(A,'double') | ndims(A)>2,
   error('Invalid argument.');
end;

B = A;

%end .. mirror
