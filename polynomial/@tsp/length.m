function l = length(T)
%LENGTH  Length of two-sided polynomial
%   
% L = LENGTH(T) returns the length of the two-sided polynomial matrix T.
% It is equivalent to MAX(SIZE(T)) for a nonempty two-sided polynomial
% matrix and 0 for an empty  matrix.

%       Author(s): S. Pejchova  21-07-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 21-Jul-1999 16:48:34   $

% Effect on other properties:
% L is a standard Matlab integer.

row = T.s(1);
col = T.s(2);
if row==0|col==0,
   l=0;
else,
   l=max(row,col);
end;

%end .. @tsp/length

