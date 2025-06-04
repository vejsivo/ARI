function Af = rot90(A,k);
%ROT90  Rotate a polynomial matrix by 90 degrees
%
% ROT90(A) is the 90 degree counter-clockwise rotation of matrix A.
% ROT90(A,K) is the K*90 degree rotation of A. K is an integer.
%
% See also POL/FLIPLR, POL/FLIPUD.

%       Author(s): M. Hromcik, M. Sebek 24-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 06-Mar-1998 09:59:59   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00  J.Jezek  $

% Effect on other properties:
% Af.u: UserData are deleted.
 
if nargin==1, k = 1;
else
   if length(k)~=1 | ~isa(k,'double') | ~isreal(k) | floor(k)~=k,
      error('Invalid 2nd argument; must be scalar integer.');
   end;
end;
 
k = rem(k,4);
if k < 0, k = k + 4; end;

Af.d = A.d;

switch k
 case 1,
   Af.s = fliplr(A.s);
   Af.c = flipdim( permute(A.c, [2,1,3]), 1 );
   
 case 2,
   Af.s = A.s;
   Af.c = flipdim( flipdim(A.c,1), 2 );
   
 case 3,
   Af.s = fliplr(A.s);
   Af.c = flipdim( permute(A.c, [2,1,3]), 2 );
 
 otherwise
   Af.s = A.s;
   Af.c = A.c;
end;	%switch
 
Af.v = A.v; 
Af.h = A.h;
Af.u = [];
Af.version = 3.0;

Af = class(Af,'pol');

%end .. @pol/rot90
