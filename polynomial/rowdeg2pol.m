function Y = rowdeg2pol(X,r)
%ROWDEG2POL Conversion of a polynomial matrix in row degree expanded form into the POL object.
%
%  Y = ROWDEG2POL(X,r) accepts a 3D-array describing a polynomial matrix in such a way 
%  that X(:,:,1) is a leading row coefficients matrix of the polynomial matrix Y, 
%  X(:,:,2) is a matrix whose rows are formed by the coefficients of the powers deg(X,'row')-1 
%  in the corresponding row of the polynomial matrix Y and so on.
%
%  The second input argument r is a vector specifiing row degrees of the resulting polynomial matrix X.
%
%  The function is inverse to POL2ROWDEG.
%
%  See also: COLDEG2POL

% Author(s):    Z. Hurak   09-24-2002
% Copyright (c) 2002 by PolyX, Ltd.
% $Revision: 1.0.0 $  $Date: 0-0-0 $

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

if nargin > 2 
    error('Too many input arguments.')
end

if nargin < 2
    error('Not enough input arguments.')
end

[xs1,xs2,xs3] = size(X);

%=============================================================================
% Reconstructing the original polynomial matrix
%=============================================================================

Y = zeros(xs1,xs2,xs3);         % prealocation for the coefficients 
                                % of the polynomial matrix                            
shiftBy = xs3-r-1;              % by how many steps must the row be shifted
for ii = 1:xs1                  % row by row
    Y(ii,:,shiftBy(ii)+1:end) = X(ii,:,1:end-shiftBy(ii));
end

%=============================================================================
% Building the polynomial matrix: 
%=============================================================================

Y = Y(:,:,end:-1:1);            
Y = pol(Y);