function B = bezout(p,q,n)
%BEZOUT   Bezoutian matrix of two polynomials
%
% Given two scalar polynomials P and Q of degree N, the instruction
%
%   B = BEZOUT(P,Q)
%
% computes the NxN symmetric Bezoutian matrix B whose entries B(i,j) satisfy
%
%  P(s)*Q(t) - Q(s)*P(t)
%  --------------------- = sum_{i=1}^N sum_{j=1}^N B(i,j)*s^(i-1)*t^(j-1)
%          s - t
%
% A third argument may be specified to force the degree of P and Q
% to a value greater than N. In this case, some leading coefficients
% in the polynomials are artificially set to zero.

% Author: Didier Henrion, January 25, 2000.
% Copyright 2000 Polyx, Ltd.
% Modified by J. Jezek 08-Aug-2001  arg check 

if nargin < 2,
 error('Not enough input arguments.');
end;

eval('p = pol(p); q = pol(q);', 'error(peel(lasterr));');

if (length(p) > 1 ) | (length(q) > 1),
  error('Scalar polynomials only.');
end;

[tv,v,p,q] = testvp(p,q);
if tv==0,
   warning('Inconsistent variables.');
elseif tv==2,
   error('Inconsistent variables.');
end;

[th,h,p,q] = testhp(p,q,v);
if ~th,
   warning('Inconsistent sampling periods.');
end;

% extract coefficient vectors
if nargin < 3,
 n = max(deg(p),deg(q));
end;

if (n < deg(p)) | (n < deg(q)),
 error('3rd argument must be greater than the degree of P and Q.');
elseif isinf(n) | (n == 0),
 B = []; return;
end;

p = p{:}; p = [p zeros(1, n+1-length(p))];
q = q{:}; q = [q zeros(1, n+1-length(q))];

% build Bezoutian entry-wise
B = zeros(n);
for i = 1:n,
 for j = i:n,
   bij = 0;
   for k = i+j-n-1:i-1,
    if (k >= 0) & (k <= n) & (i+j-k-1 >= 0) & (i+j-k-1 <= n),
     bij = bij + p(i+j-k)*q(k+1) - q(i+j-k)*p(k+1);
    end;
   end;
   if j > i,
    B(i,j) = bij; B(j,i) = bij;
   else
    B(i,i) = bij;
   end;
 end;
end;

%end .. bezout
