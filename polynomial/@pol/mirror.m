function B = mirror(A)
%MIRROR   Mirror image of polynomial
%
% B = MIRROR(A) where A is polynomial in s, p, z or z^-1,
% return the mirror image of A:
%         B(s) = A(-s)   or    B(z) = A(z^-1)
%
% See also POL/CTRANSPOSE.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.

B = A; Av = A.v;
if strcmp(Av,'z'), B.v = 'z^-1';
elseif strcmp(Av,'z^-1'), B.v = 'z';
elseif strcmp(Av,'s') | strcmp(Av,'p'),
   for i = 1:B.d+1,
      B.c(:,:,i) = (-1)^(i-1)*B.c(:,:,i);
   end;
elseif ~isempty(Av),
   error('Invalid variable symbol; must be s, p, z or z^-1.');
end;

%end .. @pol/mirror
