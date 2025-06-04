function Act = ctranspose(A)
%CTRANSPOSE   Conjugated transpose two-sided polynomial
%                  A'
% The commands
%    AT = A'
%    AT = CTRANSPOSE(A)
% return the conjugated transpose of the tsp matrix A:
%    AT(z) = A'(z^-1)
%
% See also TSP/TRANSPOSE, TSP/CONJ.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.
%        $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                          $Date: 22-May-2000  12:00:00  $

PP = ctranspose(A.p);
Act = tsp(PP); Act.h = A.h;
Act = shift(Act,-A.o);

%end .. @tsp/ctranspose
