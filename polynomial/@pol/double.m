function A = double(P)
%DOUBLE  Convert polynomial to double (standard Matlab matrix)
%             A = DOUBLE(P)
%
% The command
%     A = DOUBLE(P)
% converts a polynomial matrix P to double (standard MATLAB matrix)
% if possible (if it is empty or constant polynomial matrix).
% If not possible, error.
%
% See also POL/DECLASS.

%        Author: J. Jezek   11-10-99
%        Copyright (c) 1999 by Polyx, Ltd.

if isempty(P) | P.d<0,
   A = zeros(P.s(1),P.s(2));
elseif P.d==0,
   A = P.c(:,:,1);
else
   error('Argument is not convertible to double.');
end;

%end .. @pol/double

