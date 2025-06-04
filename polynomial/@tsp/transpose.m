function At = transpose(A)
%TRANSPOSE  (.')  Transpose two-sided polynomial
%                     A.'
%
% AT = A.'  or  AT = TRANSPOSE(A) returns the non-conjugated transpose
% of the tsp matrix A.
%
% See also TSP/CTRANSPOSE.

%       Author:  J. Jezek,  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $ $Date: 29-Sep-1999  13:00:00  $
%                        $Date: 22-May-2000  12:00:00  $

PP = transpose(A.p);
At = tsp(PP); At.h = A.h;
At = shift(At,A.o);

%end .. @tsp/transpose
