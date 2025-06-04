function y = edgetest(varargin) 
%EDGETEST  Test for robust stability of a polytope of polynomials
% 
% [Y =]EDGETEST(p0,p1,...,pk,QBOUNDS) performs a test for robust stability of 
% a polynomial family with affine dependence of the coefficientson the vector 
% of uncertain parameters. The uncertain parameters are bounded by a box 
% in the space of parameters. The polynomial family is then given by
%
% P(s,q) = P0 + q1*p1 + q2*p2 + ... + qk*pk, where
%
% q1_min <= q1 <= q1_max
% q2_min <= q2 <= q2_max
% .
% .
% .
% qk_min <= qk <= qk_max
%
% which is expressed by writing: 
%
% Qbounds = [q1_min, q1_max; q2_min, q2_max; ... ; qk_min, qk_max];
%
% The function EDGETEST returns 1 if the polynomial family is robustly stable,
% i.e. if the polynomials coresponging to the edges of the uncertain parameter
% box are stable. Otherwise, the function EDGETEST returns zero.
%
% See also: PTOPPLOT, PTOPEX

% Author(s):    Z. Hurak   11-30-2000
% Copyright (c) 2000 by PolyX, Ltd.
% $Revision: 1.0.4 $    $Date: 12-21-2000 $
%-----------------------------------------------------------------------------------------

global PGLOBAL;
eval('PGLOBAL.FORMAT;',...
   'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

%=============================================================================
%Tests for proper inputs 
% (only partial test because ptopex takes care for the rest):
%=============================================================================
ni = nargin; nr = ni - 2; y = [];

if ~ni
   error('Not enough input arguments.');  
end

if ni <= 2
   eval('y = isstable(varargin{1});',...
      'error(''The input argument is not a polynomial matrix.'');')
   return
end

%=============================================================================
%Create the 'names' for the input parameters (p0 for the nominal polynomial,
%p1,...,pk for the remainng polynomials, Qbounds for the uncertain parameter
%bounding set. And then assing the values of the input parameters to these
%names:
%=============================================================================

extremPol0 = varargin{1};
for ii = 1:nr
   eval(['extremPol',num2str(ii),'=varargin{ii+1};']);
end
Qbounds = varargin{ni};
   
%=============================================================================
%Internal call of PTOPEX function to generate the set of extreme polynomials:
%=============================================================================

extreme = [];
inpPolString = 'extremPol0,';
for ii = 1:nr
   newPolString = sprintf('extremPol%d,',ii);
   inpPolString = [inpPolString,newPolString];
end
eval(['extremePolynomials = ptopex(',inpPolString,'Qbounds);']);

%=============================================================================
%Stability test of the first extreme polynomial:
%=============================================================================

if ~isstable(extremePolynomials(1))
   y = 0;
   y = logical(y);
   return;
end

%=============================================================================
%Internal call of STABINT function to test the stability of nr*2^(nr-1) edges:
%=============================================================================

np = 2^nr;

for ii = 2:np
   [rmin, rmax] = stabint(extremePolynomials(ii-1),extremePolynomials(ii));
   if rmax ~= inf
      y = 0;
      y = logical(y);
      return
   end
end

[rmin, rmax] = stabint(extremePolynomials(np),extremePolynomials(1));
if rmax ~= inf
   y = 0;
   y = logical(y);
   return
end
   
y = 1;
y = logical(y);

% end ../edgetest

  
  
