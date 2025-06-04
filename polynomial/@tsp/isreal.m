function flag = isreal(T,tol)
%ISREAL   Test if two-sided polynomial is real
%
% ISREAL(A) returns 1 if all the entries of the tsp
% matrix A have zero imaginary part and 0 otherwise.
% 
% See also TSP/REAL, TSP/IMAG, I, J.

%    Author: J. Jezek 11-8-99
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision $  $ Date 28-Feb-2003 $

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double'),
      eror('Invalid 2nd argument.');
   end;
end;

flag = isreal(T.p);

%end .. @tsp/isreal
