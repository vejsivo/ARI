function flag = isreal(A,tol)
%ISREAL  Test if polynomial is real
%
% ISREAL(A) returns 1 if all the entries of the polynomial 
% matrix A have zero imaginary part and 0 otherwise.
% 
% See also POL/REAL, POL/IMAG, I, J.

%    Author: D. Henrion, June 19, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision $  $ Date 28-Feb-2003  J.Jezek  tol $

if nargin==2 & ~isempty(tol),
   if ~isa(tol,'double'),
      error('Invalid 2nd argument.');
   end;
end;

flag = isreal(A.c);

%end .. @pol/isreal


