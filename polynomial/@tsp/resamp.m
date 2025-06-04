function B = resamp(A,k,h)
%RESAMP    Resample two-sided polynomial
%
% Let A be a two-sided polynomial
%    A = ... + A(-n)*z^n + ... + A(0) + ... + A(n)*z^-n + ...
% whose n-th coefficient is A(n). The command
%   B = RESAMP(A,K,H)
% where scalar integers K,H are resampling period, K>=1, and
% resampling phase, H>=0, (defaults being K=2, H=0 ),
% returns two-sided polynomial B, with n-th coefficient
%  B(n) = A(K*n+H). The meaning of variable 'z' has changed:
%  new var = old var to the K-th power.
% The sampling period  B.h  of B is K-times greater than A.h .
%
% For example:
%    A = 0.03125*z^5 + 0.0625*z^4 + 0.125*z^3 + 0.25*z^2 + 0.5*z +
%        1 + 0.8*z^-1 + 0.64*z^-2 + 0.512*z^-3 + 0.4096*z^-4 +
%        0.32768*z^-5
%    RESAMP(A,3,0) = 0.125*z + 1 + 0.512*z^-1
%    RESAMP(A,3,1) = 0.03125*z^2 + 0.25*z + 0.8 + 0.4096*z^-1
%    RESAMP(A,3,2) = 0.0625*z^2 + 0.5*z + 0.64 + 0.32768*z^-1
%
% H can also be, instead of scalar, a vector. In such
% a case, the result is a cell vector.
%
% See also POL/RESAMP, RDF/RESAMP, LDF/RESAMP.

%        Author:  J. Jezek  08-Dec-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 26-May-2000 $
%                      $ Date 04-Oct-2000 $
%                      $ Date 31-Oct-2000 $
%                      $ Date 01-Feb-2001 $
%                      $ Date 28-Feb-2003 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'tsp'),
   error('Some argument but not 1st is invalidly two-sided polynomial.');
end;

if ni<2 | isempty(k),
   k = 2;
else
   if ~isa(k,'double') | length(k)~=1 | ~isreal(k) | ...
         floor(k)~=k | k<=0,
      error('Invalid resampling period.');
   end;
end;

if ni<3 | isempty(h),
   h = 0;
else
   if ~isa(h,'double'),
      error('Invalid resampling phase.');
   end;
end;

oo = floor(A.o/k); o = k*oo;
Ap = shift(A.p,A.o-o,A.v);
BB = 0;
eval('BB = resamp(Ap,k,h);','error(peel(lasterr));');

lh = length(h);
if lh==1,
   B = tsp(BB); B.h = A.h;
   if ~isempty(B.h), B.h = B.h*k;
   end;
   B = shift(B,oo);
else
   B = cell(size(h));
   for i = 1:lh,
      Q = tsp(BB{i}); Q.h = A.h;
      if ~isempty(Q.h), Q.h = Q.h*k;
      end;
      B{i} = shift(Q,oo);
   end;
end;

%end .. @tsp/resamp

