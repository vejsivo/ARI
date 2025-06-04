function [Ai,Ar] = imag(A,tol)
%IMAG    Imaginary part of two-sided polynomial
%             AI = IMAG(A)
%
% IMAG(A)  is the imaginary part of A.
% [AI,AR] = IMAG(A) returns also the real part.
%
% See also TSP/REAL.

%       Author: J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 22-May-2000  12:00:00  $
%                         $Date: 28-Feb-2003  tol       $

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double'),
      error('Invalid 2nd argument.');
   end;
end;

PP = imag(A.p);
Ai = tsp(PP); Ai.h = A.h;
Ai = shift(Ai,A.o);

if nargout==2,
   PP = real(A.p);
   Ar = tsp(PP); Ar.h = A.h;
   Ar = shift(Ar,A.o);
end;

%end .. @tsp/imag
