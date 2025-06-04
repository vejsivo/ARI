function y = isrowred(M)
%ISROWRED Test for row reducedness of a polynomial matrix.
%
%[Y=]ISROWRED(M) returns 1 if the polynomial matrix is row reduced 
%            and returns 0 otherwise.
%
% See also: LCOEF, ROWRED.

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
% Tests for row reducedness: 
%=============================================================================

[m,n] = size(M);
r = rank(lcoef(M,'row'));

y = r==min(m,n);

