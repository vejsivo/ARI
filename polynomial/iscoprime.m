function r = iscoprime(P,tol)
%ISCOPRIME   Test if constant matrix is coprime
%
% The command
%    R = ISCOPRIME(F[,TOL])
% for constant matrix F returns always 1.
%
% This macro exists only for completeness.

%      Author:  J. Jezek, 02-Feb-2003
%      Copyright(c) 2003 by Polyx, Ltd.

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(P,'double') | ndims(P)>2,
   error('Invalid 1st argument.');
end;

r = logical(1);

%end .. iscoprime
