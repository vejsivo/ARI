function [y,emax] = isdiscstable(x,c)
%ISDISCSTABLE Test for robust stability of a disc polynomial.
%
% Y = ISDISCSTABLE(X,C) returns 1 (true) if the disc polynomial given 
% by the nominal polynomial X and the vector of radii of uncertainty C 
% is robustly stable. It returns 0 (false) if the disc polynomial is not 
% robustly stable. Both continuous and discrete stability is tested.
%
% [Y,EMAX] = ISDISCSTABLE(X,C) returns an additional output EMAX such 
% that the disc polynomial with uncertainty radii less than EMAX*C is 
% robustly stable.
%
% See also: KHARIT.

% Author(s):    Z. Hurak   06-22-2001
% Copyright (c) 2001 by PolyX, Ltd.
% $Revision: 1.0.0 $  $Date: 06-22-2001 $

%-----------------------------------------------------------------------------

global PGLOBAL;
eval('PGLOBAL.FORMAT;',...
   'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

%=============================================================================
% Tests for proper inputs: 
%=============================================================================

if nargin ~= 2
    error('Two input arguments required.')
end

if ~isa(x,'pol')      % dead branch for polynomial matrices...
    error('The first input argument should be a polynomial.')
end

if ~isnumeric(c)
    error('The second input argument should be a numeric array.')
end

[mx,nx] = size(x);
degx = deg(x);
leC = length(c);

if ~(mx == 1 & nx == 1) 
    error('Current version supports scalar polynomials only.')
end

if degx+1 ~= leC
    error('Inconsistend degree of the polynomial and the dim of the vector of radii.')
end

%=============================================================================
% Tests for proper outputs: 
%=============================================================================

if nargout > 2
    error('Two many output arguments.')
end

%=============================================================================
% Tests for Hurwitzness or Schurness of the nominal polynomial: 
%=============================================================================

if ~isstable(x)
    y = 0;
    return
end

%=============================================================================
% Computation of the Hinf norms of the two rational functions (continuous case): 
%=============================================================================
if isct(x)    
	qSign = ceil(leC/4);
	qImag = ceil(leC/2);
	signPattern = repmat([-1 -1 1 1],1,qSign);
	imagPattern = repmat([1 j],1,qImag);
    imagPattern = imagPattern(1:leC);
	
	sign1 = [1 signPattern]; sign1 = sign1(1:leC);
	sign2 = [1 1 signPattern]; sign2 = sign2(1:leC);
	
	cImag = c.*imagPattern;
	gamma1coef = cImag.*sign1;
	gamma2coef = cImag.*sign2;
	
	gamma1 = pol(gamma1coef, degx);
	gamma2 = pol(gamma2coef, degx);
	
	g1 = hinfnorm(gamma1,x);
	g2 = hinfnorm(gamma2,x);
	
	y = (g1<1) & (g2<1);
    
	if nargin == 2
        emax = min(1/g1,1/g2);
	end
end


%=============================================================================
% Computation of the infimum of the polynomial on unit circle (discrete case): 
%=============================================================================

if isdt(x)
    N = 4096;           % number of Fourier points
    xx = x{:};
    XX = abs(fft(xx,N));
    minXX = min(XX);
    sumC = sum(c);
    if nargout == 2
        emax = minXX/sumC;
    end
    y = sumC<minXX;
    
end

% end of ..\isdiscstable.m