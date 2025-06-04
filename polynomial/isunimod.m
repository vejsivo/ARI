function r = isunimod(P,arg2)
%ISUNIMOD  Test is constant matrix is unimodular
%
% The command
%    R = ISUNIMOD(P[,TOL])
% returns 1 if the square constant matrix is nonsingular
% and 0 otherwise.
%
% This macro exists only for completeness.
% See also POL/ISUNIMOD, TSP/ISUNIMOD.

%      Author:  J. Jezek, 19-Jul-2001
%      Copyright(c) 2001 by Polyx, Ltd.
%      $ Revision $  $ Date 25-Jul-2002 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(P,'double') | ndims(P)>2,
   error('Invalid 1st argument.');
end;

if ni==1,
   eval('r = ~issingular(P);', 'error(peel(lasterr));');
else
   eval('r = ~issingular(P,arg2);', 'error(peel(lasterr));');
end

%end .. isunimod
