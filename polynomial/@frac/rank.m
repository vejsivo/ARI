function r = rank(A,arg2,arg3)
%RANK  Rank of fraction
%
% The command
%    RANK(A) 
% provides an estimate of the number of linearly independent rows 
% or columns of the fraction A (left-den fraction, right-den
% fraction or scalar-den fraction).
%
% An optional argument MET may specify the numerical method to be used,
% like for polynomials.
%
% An optional argument TOL may specify the zeroing tolerance to be used
% instead of the standard one.
%
% See also  POL/RANK.

%       Author: J. Jezek  02-Feb-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 14-Oct-2002 $

na = nargin;
if na==1,
   r = rank(A.num); return;
end;

if ~isa(arg2,'double') & ~isa(arg2,'char'),
   error('Invalid 2nd argument.');
end;
if na==3,
   if ~isa(arg3,'double') & ~isa(arg3,'char'),
      error('Invalid 3rd argument.');
   end;
end;

if na==2,
   eval('r = rank(A.num,arg2);','error(peel(lasterr));');
else
   eval('r = rank(A.num,arg2,arg3);','error(peel(lasterr));');
end;

%end .. @frac/rank

