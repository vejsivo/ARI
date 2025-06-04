function [s1,s2,s3] = size(R,d)
%SIZE  Size of a fraction
%
% Given an M-by-N fraction R, the command
%    D = SIZE(R)
% returns the two-element row vector D = [M N] containing the 
% numbers of rows and columns of the matrix.
% The command
%    [M,N] = SIZE(R)
% returns the numbers of rows and columns in separate output variables.
% The command
%    [M,N,D] = SIZE(R)
% returns the number of rows N, the number of columns M and
% the degree of fraction R.
%
% M = SIZE(R,1) returns just the number of rows.
% N = SIZE(R,2) returns just the number of columns.
% D = SIZE(R,3) returns just the degree of fraction R.

%       Author:  J. Jezek  24-Jan-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Feb-2003 $

ni = nargin;  no = nargout;
if ni==2 & (~isa(d,'double') | length(d)~=1 | ...
      ~isreal(d) | round(d)~=d),
   error('Invalid 2nd argument.');
end;
row = R.s(1);
col = R.s(2);
if ni==1,
   s1 = row; s2 = col;
   if no<=1,
      s1 = [row,col];
   elseif no==3,
      s3 = deg(R);
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
      s1 = deg(R);
   otherwise,
      error('Invalid 2nd argument; must be 1, 2 or 3.');
   end;
end;

%end .. @frac/size

