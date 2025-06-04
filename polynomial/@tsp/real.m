function [Ar,Ai] = real(A,tol)
%REAL    Real part of two-sided polynomial
%             AR = REAL(A)
%
% REAL(A) is the real part of A.
% [AR,AI] = REAL(A) returns also the imag part.
%
% See also TSP/IMAG.

%     Author: J. Jezek  11-8-99
%     Copyright (c) 1999 by Polyx, Ltd.
%     $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                       $Date: 22-May-2000  12:00:00  $
%                       $Date: 28-Feb-2003   tol      $

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double'),
      error('Invalid 2nd argument.');
   end;
end;

PP = real(A.p);
Ar = tsp(PP); Ar.h = A.h;
Ar = shift(Ar,A.o);

if nargout==2,
   PP = imag(A.p);
   Ai = tsp(PP); Ai.h = A.h;
   Ai = shift(Ai,A.o);
end;

%end .. @tsp/real
