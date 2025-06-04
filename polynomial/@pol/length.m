function l = length(P)
%LENGTH  Length of polynomial
%   
% L = LENGTH(P) returns the length of the polynomial matrix P. It is
% equivalent to MAX(SIZE(P)) for a nonempty polynomial matrix and
% 0 for an empty polynomial matrix.

%       Author(s): S. Pejchova, M. Sebek 08-09-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 8-Sep-1998 10:22:34   $

% Effect on other properties:
% L is a standard Matlab integer.

narginchk(1,1);
% error(nargchk(1,1,nargin));	%REMOVED IN NEW MATLABS
if nargout>1, error('Too many output arguments.'); end;

row = P.s(1);
col = P.s(2);
if row==0|col==0,
   l=0;
else,
   l=max(row,col);
end;

%end .. @pol/length
