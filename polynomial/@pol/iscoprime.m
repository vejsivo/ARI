function r = iscoprime(P,tol)
%ISCOPRIME   Test if polynomial is coprime
%
% The command
%    R = ISCOPRIME(F[,TOL])
% for polynomial F returns always 1.
%
% This macro exists only for completeness.

%      Author:  J. Jezek, 02-Feb-2003
%      Copyright(c) 2003 by Polyx, Ltd.

if ~isa(P,'pol'),
   error('Invalid 1st argument.');
end;

r = logical(1);

%end .. @pol/iscoprime
