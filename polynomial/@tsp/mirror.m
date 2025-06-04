function B = mirror(A)
%MIRROR   Mirror image of two-sided polynomial
%
% The command
%      B = MIRROR(A)
% returns the mirror image of the tsp matrix A:
%              B(z) = A(z^-1)
%
% See also TSP/CTRANSPOSE.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) by Polyx, Ltd.
%        $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                          $Date: 13-Oct-1999  12:00:00  $

B = A;
Bp = B.p; Br = B.r;
if isempty(Bp) | Br<0, return;
end;
Bpc = Bp.c; Bpc = flipdim(Bpc,3);
B.p = pol(Bpc(:,:),Br); B.p.v = 'z';
B.d = -A.t; B.t = -A.d;
B.o = -A.o-Br;

%end .. @tsp/mirror

