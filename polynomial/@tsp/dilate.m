function B = dilate(A,k,h,var)
%DILATE     Dilate two-sided polynomial
%
% Let A be a two-sided polynomial
%    A(z^-1) = ... + A(-n)*z^n + ... + A(0) + ... + A(n)*z^-n + ...
% whose n-th coefficient is A(n). The command
%    B = DILATE(A,K,H)
% where scalar integers K,H are dilating period, K>=1, and dilating
% phase, H>=0, returns two-sided polynomial
%   B(z^-1) = z^-H * A(z^-K). Note that only such coefficients of B
% are nonzero which can be expressed B(K*n+H) = A(n) with some n.
% The arguments K,H are optional, the defaults being  K=2, H=0.
%
% The meaning of the variable 'z' has changed:
%  old var = new var to the K-th power
% The sampling period of B is K-times less than that of A.
%
% See also POL/DILATE, TSP/RESAMP.

%        Author:  J. Jezek 08-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 26-May-2000 $
%                      $ Date 04-Oct-2000 $
%                      $ Date 28-Feb-2003 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
eval('A = tsp(A);', ...
   'error(''Invalid 1st argument.'');');

if ni<2 | isempty(k),
   k = 2;
else
   if ~isa(k,'double') | length(k)~=1 | ~isreal(k) | floor(k)~=k ...
         | k<1,
      error('Invalid dilating period.');
   end;
end;

if ni<3 | isempty(h),
   h = 0;
else
   if ~isa(h,'double') | length(h)~=1 | ~isreal(h) | floor(h)~=h ...
         | h<0,
      error('Invalid dilating phase.');
   end;
end;

PP = dilate(A.p,k);
B = tsp(PP); B.h = A.h;
if ~isempty(B.h), B.h = B.h/k;
end;
B = shift(B,k*A.o-h);

%end .. @tsp/dilate
