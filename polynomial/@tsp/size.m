function [s1,s2,s3] = size(T,d)
%SIZE  Size of a two-sided polynomial matrix
%
% Given an M-by-N two-sided polynomial matrix T, the command
%    D = SIZE(T)
% returns the two-element row vector D = [M N] containing the 
% numbers of rows and columns of the matrix.
% The command
%    [M,N] = SIZE(T)
% returns the numbers of rows and columns in separate output variables.
% The command
%    [M,N,Dg] = SIZE(T)
% returns the number of rows, of columns and the degree of T.
% 
% M  = SIZE(T,1) returns just the number of rows.
% N  = SIZE(T,2) returns just the number of columns.
% Dg = SIZE(T,3) returns just the degree of T.

%       Author(s): S. Pejchova  21-07-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 20-Sep-1999 11:40:34   $
%                         $Date: 28-Feb-2003  J.Jezek $

% Effect on other properties:
% [s1,s2,s3] is a standard Matlab vector.

ni = nargin;  no = nargout;
if ni==2 & isa(d,'tsp'),
   error('Invalid 2nd argument, must be 1, 2 or 3.');
end;
row = T.s(1);
col = T.s(2);
if ni==1,
   s1 = row; s2 = col;
   if no<2,
      s1 = [row,col];
   elseif no==3,
      s3 = deg(T);
   end;
elseif no > 1,
   error('Too many output arguments.');
else,
   switch d,
   case 1,
      s1 = row;
   case 2,
      s1 = col;
   case 3,
      s1 = deg(T);
   otherwise,
      error('Invalid 2nd argument; must be 1, 2 or 3.');
   end;
end;

%end .. @tsp/size
