function r = rank(A,arg2,arg3)
%RANK  Rank of matrix-den fraction
%
% The command
%    RANK(A) 
% provides an estimate of the number of linearly independent rows 
% or columns of the matrix-den fraction A.
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
%       $ Revision $  $ Date 26-Apr-2000 $

eval('A = mdf(A);','error(peel(lasterr));');
met = ''; tol = [];
na = nargin;
if na>=2,
   if isa(arg2,'double') & ...
         (isempty(arg2) | (length(arg2)==1 & arg2>=0)),
      tol = arg2;
   elseif isa(arg2,'char') & ...
         (isempty(arg2) | strcmp(arg2,'fft') | strcmp(arg2,'syl')),
      met = arg2;
   else
      error('Invalid 2nd argument.');
   end;
end;
if na==3,
   if isa(arg3,'double') & ...
         (isempty(arg3) | (length(arg3)==1 & isreal(arg3) & arg3>0)),
      tol = arg3;
   elseif isa(arg3,'char') & ...
         (isempty(arg3) | strcmp(arg3,'fft') | strcmp(arg3,'syl')),
      met = arg3;
   else
      error('Invalid 3rd argument.');
   end;
end;

if A.frac.s(1)>=A.frac.s(2), F = ldf(A,tol);
else F = rdf(A,tol);
end;
r = rank(F.n,met,tol);

%end .. @mdf/rank
