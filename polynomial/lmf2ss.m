function [a,b,c,d] = lmf2ss(N,D,tol)
%LMF2SS  Converts a left matrix fraction to an observer-form realization {A,B,C,D}
%
% Given two polynomial matrices N and D such that D is column reduced
% the command
%    [a,b,c,d] = LMF2SS(N, D [,TOL])
% returns the (generalized) observer-form realization (a,b,c,d) of the
% system with transfer matrix H, that is,
%     H(x) = D^-1(x) N(x) = c (xI-a)^-1 b + d(x),
% where x is 's' (continuous-time) or 'z' (discrete-time).
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: RMF2SS.

%    Author: D. Henrion, M. Sebek, R.C.W. Strijbos, November 13, 1998.
%    Copyright 1998 by Polyx, Ltd.

if nargin < 2,
   error('Not enough input arguments.');
end;
if nargin==3 & ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

eval('N = pol(N); D = pol(D);', ...
   'error(peel(lasterr));');

% Call RMF2SS.

if nargin == 2
   eval('[a,b2,c2,d] = rmf2ss(N.'', D.'');', ...
      'error(peel(lasterr));');
else
   eval('[a,b2,c2,d] = rmf2ss(N.'', D.'', tol);', ...
      'error(peel(lasterr));');
end
a = a'; b = c2'; c = b2'; d = d.';

%end .. lmf2ss

