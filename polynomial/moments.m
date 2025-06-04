function ms = moments(p)
%MOMENTS Moment of roots of a polynomial.
%
% MOMENTS(p) returns the first n moments of roots s(i) of an n-th degree polynomial 
% p(v)=a(n)*v^n+a(n-1)*v^(n-1)+...+a(1)*v+a(0), where the moments of the roots are defined as
%
% m(1) = s(1)   + s(2)   + ... + s(n),
% m(2) = s(1)^2 + s(2)^2 + ... + s(n)^2,
% m(3) = s(1)^3 + s(2)^3 + ... + s(n)^3,
% .
% .
% .
% m(n) = s(1)^n + s(2)^n + ... + s(n)^n,
%
% The first n moments of an n-th degree monic polynomial uniquely determine its coefficients.
%
% See also: ROOTS

% Author(s):    Z. Hurak   11-27-2001
% Copyright (c) 2001 by PolyX, Ltd.
% $Revision: 1.0.0 $    $Date: 11-27-2001 $
%-----------------------------------------------------------------------------------------

global PGLOBAL;
eval('PGLOBAL.FORMAT;',...
   'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

%====================================================================
%Tests for proper inputs:
%====================================================================

if nargin ~= 1
    error('One input argument required.')
end

if ~isa(p,'pol')
    error('Input argument should be a POL object.')
end

[nr,nc] = size(p);

if nc*nr ~= 1
    error('Input argument should be a scalar polynomial only.')
end

%====================================================================
% Building and solving a linear system:
%====================================================================

dp = deg(p);                               % degree of a polynomial p
pCoeffs = fliplr(p{:});                    % reverse order of coeffs.
pCoeffs = pCoeffs/pCoeffs(1);              % make the polynomial monic
A = toeplitz(pCoeffs(1:end-1),eye(1,dp));  % toeplitz matrix
b = -(1:dp).*pCoeffs(2:end);               % right hand side for Ax=b
ms = A\b';

% end of ..\moments.m