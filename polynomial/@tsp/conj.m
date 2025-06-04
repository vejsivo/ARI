function Ac = conj(A)
%CONJ    Complex conjugate of two-sided polynomial
%             Ac = conj(A)
%
% The commmand
%    CONJ(A) 
% returns the complex conjugate of A, that is
%    CONJ(A) = REAL(A) - i*IMAG(A)
%
%   See also TSP/CTRANSPOSE, TSP/REAL, TSP/IMAG, I, J.

%     Author: J. Jezek  11-8-99
%     Copyright (c) 1999 by Polyx, Ltd.
%     $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                       $Date: 22-May-2000  12:00:00  $

PP = conj(A.p);
Ac = tsp(PP); Ac.h = A.h;
Ac = shift(Ac,A.o);

%end .. @tsp/conj
