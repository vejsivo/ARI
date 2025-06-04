function y = iscolred(M)
%ISCOLRED Test for column reducedness of a polynomial matrix.
%
%[Y=]ISCOLRED(M) returns 1 if the polynomial matrix is column reduced 
%            and returns 0 otherwise.
%
% See also: LCOEF, COLRED.

% Author(s):    Z. Hurak   06-20-2001
% Copyright (c) 2001 by PolyX, Ltd.
% $Revision: 1.0.0 $  $Date: 06-20-2001 $

%---------------------------------------------------------------------------------

global PGLOBAL;
eval('PGLOBAL.FORMAT;',...
   'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

%=============================================================================
% Tests for proper inputs: 
%=============================================================================

% no tests, because LCOEF works with constant matrices too.

%=============================================================================
% Tests for column reducedness: 
%=============================================================================

[m,n] = size(M);
r = rank(lcoef(M,'col'));

y = r==min(m,n);

