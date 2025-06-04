function S = jury(p, d, flag)
%JURY  Create the Jury matrix corresponding to a polynomial
%
% The command
%    J = JURY(P,K)
% creates the constant Jury matrix J of dimension (K-1)-by-(K-1)
% corresponding to the polynomial P. If
%    P = P{0} + P{1}*v + P{2}*v^2 + ... + P{D}*v^D
% and K <= D, then
%
%    J = [ P{K} P{K-1} .. P{2} ]   [    0   0  ..   P{0} ]
%        [    0   P{K} .. P{3} ]   [           ..        ]
%        [             ..      ] - [    0 P{0} .. P{K-3} ]
%        [    0     0  .. P{K} ]   [ P{0} P{1} .. P{K-2} ]
%
% The default value of K is D.
%
% With the syntax
%     J = JURY(P,K,'rev')
% coefficients p{0}, p{1}, .., p{K} are reversed.

% Author: D. Henrion, April 10, 2000.
% Updated to 3.0 by D. Henrion, August 11, 2000.
% Copyright 2000 by Polyx, Ltd.

if nargin < 1,
   error('Not enough input arguments.');
end;
eval('p = pol(p);', 'error(peel(lasterr));');
flag = (nargin > 2);
if nargin < 2 | isempty(d),
   d = deg(p);
else
   if ~isa(d,'double') | length(d)~=1 | ~isreal(d),
      error('Invalid 2nd argument.');
   end;
   if d > deg(p),
      d = deg(p);
   end;
end;

q = p{:};
if length(q) < d+1,
 q = [q zeros(1, d+1-length(q))];
end;

if flag,
 S = fliplr(hankel(q(d+1:-1:3))) - rot90(hankel(q(1:d-1)), 2);
else
 S = fliplr(hankel(q(3:d+1))) - rot90(hankel(q(d-1:-1:1)), 2);
end;

%end .. jury
