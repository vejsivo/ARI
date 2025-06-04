function H = hermfuji(p)
% The instruction
%
%   H = HERMFUJI(P)
%
% returns the Hermite (or Hermite-Fujiwara) matrix of real polynomial P.
% H is a symmetric matrix which is bilinear in coefficients of P.
% H is positive definite if and only if P is a stable polynomial.
%
% See also: HERMCOEF

% Written by D. Henrion, November 18, 2002
% Copyright 2002 by PolyX, Ltd.

global PGLOBAL;

eval('PGLOBAL.FORMAT;',...
     'error(''Use PINIT to initialize the Polynomial Toolbox.'');');
  
p = pol(p);
d = deg(p);
if d < 1,
 error('Input argument must be a non-constant polynomial');
end;

% stability region
switch symbol(p),
 case {'s','p'}
   S = [0 1;1 0]; 
 case {'z^-1','d'}
   S = [1 0;0 -1];
 case {'z','q'}
   S = [-1 0;0 1];
 otherwise
   error('Invalid input polynomial');
end;

% build vector [p0*p0;p1*p0;p1*p1;p2*p0;p2*p1;p2*p2..]
pcoef = zeros((d+1)*(d+2)/2,1); k = 1;
for i = 0:d
 for j = 0:i
  pcoef(k) = p{i}*p{j};
  k = k+1;
 end;
end;

% build Hermite vector [H11;H12;H22;H13;H23;H33..]
C = hermcoef(S,d);
h = C*pcoef; k = 1;
H = zeros(d);
for i = 1:d
 for j = 1:i
  H(i,j) = h(k);
  H(j,i) = h(k);
  k = k+1;
 end;
end;

