function r = ismonomod(P,arg2)
%ISMONOMOD  Test is constant matrix is monomodular
%
% The command
%    R = ISMONOMOD(P[,TOL])
% returns 1 if the square constant matrix P is nonsingular
% and 0 otherwise.
%
% This macro exists only for completeness.
% See also POL/ISMONOMOD, TSP/ISMONOMOD.

%      Author:  J. Jezek, 30-Sep-2002
%      Copyright(c) 2002 by Polyx, Ltd.

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

%end .. ismonomod
