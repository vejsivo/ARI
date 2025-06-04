function B = rot90(A,k)
%ROT90     Rotate the tsp matrix by 90 degrees
%           B = ROT90(A),    B = ROT90(A,k)
%
% ROT90(A) is the 90 degree counter-clockwise rotation of matrix A.
% ROT90(A,K) is the K*90 degree rotation of A. K is an integer.
%
% See also TSP/FLIPLR, TSP/FLIPUD.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0  $  $Date: 29-Sep-1999  13:00:00  $
%                          $Date: 24-May-2000  12:00:00  $
%                          $Date: 31-Oct-2000  13:00:00  $

if nargin==1, k=1;
else
   if ~isa(k,'double'),
      error('Invalid 2nd argument.');
   end;
end;
PP = 0;
eval('PP = rot90(A.p,k);','error(peel(lasterr));');
B = tsp(PP); B.h = A.h;
B = shift(B,A.o);

%end .. @tsp/rot90
