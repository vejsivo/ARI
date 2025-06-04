function [Y,rowdegX] = pol2rowdeg(X)
%POL2ROWDEG Decomposition of a polynomial matrix acording to row degrees.
%
%  Y = POL2ROWDEG(X) accepts a polynomial matrix X as an input and returns
%  a three-dimensional array Y such that Y(:,:,1) is a leading row
%  coefficients matrix, Y(:,:,2) is a matrix whose rows are formed 
%  by the coefficients of the powers deg(X,'row')-1 in the corresponding 
%  row of X and so on.
%
%  [Y,r] = POL2ROWDEG(X) returns also a column vector r of row degrees. 
%  This vector is necessary to reconstruct the polynomial matrix X from Y. 
%
%  See also: POL2COLDEG, ROWDEG2POL, COLDEG2POL.

% Author(s):    Z. Hurak   09-23-2002
% Copyright (c) 2002 by PolyX, Ltd.
% $Revision: 1.0.0 $  $Date: 09-23-2002 $

%-----------------------------------------------------------------------------

global PGLOBAL;
try,
    PGLOBAL.FORMAT;
catch,
    painit
end

%=============================================================================
% Tests for proper inputs: 
%=============================================================================

if nargin > 1 
    error('Too many input arguments.')
end

if nargin < 1
    error('Not enough input arguments.')
end

[xs1,xs2] = size(X);
xs3 = X.d+1;                    % "height" of the polynomial matrix X

%=============================================================================
% Building the new coefficient matrices: 
%=============================================================================

rowdegX = deg(X,'row');         % row degrees of the polynomial matrix X
rowdegX(isinf(rowdegX)) = 0;
Y = zeros(xs1,xs2,xs3);         % pre-allocation for the resulting matrix
Xc = X.c;

for ii = 1:xs1                  % row by row
    Y(ii,:,1:rowdegX(ii)+1) = Xc(ii,:,rowdegX(ii)+1:-1:1);
end

%=============================================================================
% end of ..\pol2rowdeg.m
%=============================================================================