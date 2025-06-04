function B = flipud(A)
%FLIPUD    Flip the rows of a tsp matrix in the up-down direction
%                       B = flipud(A)
%
% The command
%    FLIPUD(X) 
% returns X with its columns preserved and its rows flipped
% in the up/down direction.
%
% See also TSP/FLIPLR, TSP/ROT90.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 22-May-2000  12:00:00  $

PP = flipud(A.p);
B = tsp(PP); B.h = A.h;
B = shift(B,A.o);

%end .. @tsp/flipud
