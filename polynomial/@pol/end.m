function y = end(a,k,n)
%END Overloaded END method for a polynomial matrix
%
% If P is a polynomial matrix P = P0 + P1*s + P2*s^2 + .. + Pd*s^d
% then P(end,end) returns a scalar polynomial that is found at the last 
% last row and the last column of the polynomial matrix.
%
% See also: POL/SUBSREF, POL/SUBSASGN.

% Author(s):    Z. Hurak   07-12-2001
% Copyright (c) 2001 by PolyX, Ltd.
% $Revision: 1.0.0 $  $Date: 07-12-2001 $

%-----------------------------------------------------------------------------

if n > 2
   error('One or two indices for a polynomial matrix are allowed.');
end

sizeOfA = size(a);
y = sizeOfA(k);
if n == 1
    y = prod(sizeOfA);
end


