function B = fliplr(A)
%FLIPLR    Flip the columns of a tsp matrix in the left-right direction
%                       B = FLIPR(A)
%
% The command
%    FLIPLR(A) 
% returns A with its rows preserved and its columns flipped in the 
% left/right direction.              
%
% See also TSP/FLIPUD, TSP/ROT90.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 22-May-2000  12:00:00  $

PP = fliplr(A.p);
B = tsp(PP); B.h = A.h;
B = shift(B,A.o);

%end .. @tsp/fliplr
