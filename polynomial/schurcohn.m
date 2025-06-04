function H = schurcohn(p,n)
%SCHURCOHN    Schur-Cohn matrix of a polynomial
%
% The instruction
%   H = SCHURCOHN(P)
% computes the symmetric Schur-Cohn matrix of discrete-time scalar
% polynomial P. Matrix H is positive definite if and only if P has all
% its zeros within the unit disk.
%
% A second argument may be specified to force the degree of P
% to a value greater than N. In this case, some leading coefficients
% in the polynomial are artificially set to zero.
%
% See also: HERMFUJI, BEZOUT.

%    Author: Didier Henrion, March 8, 2000.
%    Updated to 3.0 by D. Henrion, September 1, 2000.
%    Modified by J. Jezek, Aug 19, 2001, arg checking
%    Copyright 2000 by Polyx, Ltd.

if nargin < 1,
   error('Not enough input arguments.');
end;

eval('p = pol(p);', 'error(peel(lasterr));')
if length(p) ~= 1,
   error('Scalar polynomial only.');
end; 
sp = symbol(p);
if ~(isempty(sp) | strcmp(sp,'z') | strcmp(sp,'z^-1') | ...
      strcmp(sp,'d') | strcmp(sp,'q')),
   error('Discrete-time polynomial only.');
end;

if nargin == 1
   n = deg(p);
else
   if ~isnumeric(n) | length(n)~=1 | ~isreal(n),
      error('Invalid 2nd argument.');
   end;
   if n < deg(p),
      error('Invalid 2nd argument; must be greater than deg(P).');
   end;
end;

if isinf(n) | (n == 0),
  H = []; return;
end;

r = zeros(1, n+1);
r(1:deg(p)+1) = p{:};

H = zeros(n);
for i = 0:n-1,
 for j = 0:n-1,
  for k = 0:min(i,j),
    H(i+1,j+1) = H(i+1,j+1) + p{n-i+k}*p{n-j+k} - p{i-k}*p{j-k};
  end;
 end;
end;

%end .. schurcohn
