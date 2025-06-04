function r = rank(A,arg2,arg3)
%RANK  Rank of two-sided polynomial
%
% The commands
%    RANK(A) 
%    RANK(A,'FFT') 
% provide an estimate of the number of linearly independent rows 
% or columns of the polynomial matrix A by computing the rank of A 
% evaluated at a suitable number of Fourier points.
%
% The number
%    RANK(A, TOL) or RANK(A, 'FFT', TOL) 
% is the minimal number of singular values of A evaluated at a suitable
% number of Fourier points that are larger than TOL. The commands RANK(A) 
% and RANK(A,'FFT') use the default tolerance TOL = MAX(SIZE(A))*NORM(A)*EPS).
%
% The command
%    RANK(A,'SYL') or RANK(A,'SYL',TOL) 
% estimates the rank of the polynomial matrix A by computing ranks of 
% appropriate Sylvester matrices. The algorithm is based on evaluation 
% of the rank of a constant matrices and therefore the optional input 
% argument TOL has the same meaning as for the standard MATLAB function 
% RANK.

%       Author: J. Jezek  11-8-99
%       Copyright (c) 1998 by Polyx, Ltd.

eval('A = tsp(A);','error(peel(lasterr));');

na = nargin;
if na==1,
   r = rank(A.p); return;
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
   eval('r = rank(A.p,arg2);','error(peel(lasterr));');
else
   eval('r = rank(A.p,arg2,arg3);','error(peel(lasterr));');
end;

%end .. @tsp/rank

