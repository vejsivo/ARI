function [s1,s2,s3] = size(P,d)
%SIZE  Size of a polynomial matrix
%
% Given an M-by-N polynomial matrix P, the command
%    D = SIZE(P)
% returns the two-element row vector D = [M N] containing the 
% numbers of rows and columns of the matrix.
% The command
%    [M,N] = SIZE(P)
% returns the numbers of rows and columns in separate output variables.
% The command
%    [M,N D] = SIZE(P) 
% returns the number of rows M, the number of columns N , and the degree 
% D of the polynomial matrix P.
% 
% M = SIZE(P,1) returns just the number of rows.
% N = SIZE(P,2) returns just the number of columns.
% D = SIZE(P,3) returns just the degree of polynomial matrix.

%       Author(s): S. Pejchova, M. Sebek 24-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 24-Feb-1998 18:00:34   $

% Effect on other properties:
% [s1,s2,s3] is a standard Matlab vector.

ni = nargin;
no = nargout;
narginchk(1,2);
% error(nargchk(1,2,ni));	%REMOVED IN NEW MATLABS

P = pol(P);
row = P.s(1);
col = P.s(2);
s3 = P.d;

if ni==2,
   nargoutchk(0,1);
   % error(nargchk(0,1,no));	%REMOVED IN NEW MATLABS
   if d==1,
      s1 = row;
   elseif d==2,
      s1 = col;
   elseif d==3,
      s1 = s3;
   else,
      error('Invalid 2nd argument; must be 1, 2 or 3.');
   end;
elseif no==0 | no==1,
   s1 = [row,col];
else
   s1 = row;
   s2 = col; 
end

%end .. @pol/size
