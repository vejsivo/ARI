function [c,s] = cgivens1(x,y,tol)
%CGIVENS1  Givens rotation
%
% Given two complex numbers x and y the function
%    [c,s] = cgivens1(x,y)
% computes a real number c and a complex number s such that
%    [  c    s ]
%    [ -s'   c ]
% is unitary and
%    [  c    s ] [ x ]  =  [ z ]
%    [ -s'   c ] [ y ]     [ 0 ]
% In the form
%    [c,s] = cgivens1(x,y,tol)
% x and y are set equal to zero if their magnitude
% is less than the tolerance tol. The default
% tolerance is 0.

% Author: Huibert Kwakernaak, 1997
% Copyright 1998 by PolyX Ltd.

if nargin < 2,
   error('Not enough input arguments.');
end;
if ~isnumeric(x) | length(x) ~= 1,
   error('Invalid 1st argument; must be scalar number.');
end;
if ~isnumeric(y) | length(y) ~= 1,
   error('Invalid 2nd argument; must be scalar number.');
end;
if nargin == 2, tol = 0; end
if abs(x) < tol, x = 0; end
if abs(y) < tol, y = 0; end

if x == 0
   if y == 0
      c = 1; s = 0;
   else
      c = 0; s = 1;
   end
else
   c = abs(x)/sqrt(abs(x'*x+y'*y));
   s = c*y'/x';
end

%end .. cgivens1
