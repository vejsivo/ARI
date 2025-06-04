function [num,den] = lmf2rat(N,D,tol)
%LMF2RAT  Converts a left matrix fraction to rational form
%
% Given two polynomial matrices N and D, where D is square and
% nonsingular with the same number of columns as N, the function
%    [NUM,DEN] = LMF2RAT(N, D [,TOL])
% returns the transfer function H = D^{-1} * N as polynomial
% matrices NUM and DEN that contain the numerators and the
% denominators of the entries of the transfer function matrix H, 
% respectively.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: RMF2RAT, RAT2LMF.

%    Author: R.C.W. Strijbos, November 13, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    Modified by J. Jezek, 22-Aug-2001, arg checking

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

switch nargin
case 2
   tol = PGLOBAL.ZEROING;
case 3
   if ~isempty(tol)
      if ~isa(tol,'double') | length(tol) > 1 | ...
            ~isreal(tol) | tol<0 | tol>1
         error('Invalid tolerance.');
      end
   end
otherwise
   error('Not enough input arguments.');
end
 
 eval('N = pol(N); D = pol(D);', ...
    'error(peel(lasterr));');
 eval('[N,D] = testdnd(N,D);', ...
    'error(peel(lasterr));');
 Var = ''; h = 0;
 eval('[Var,h,N,D] = testvhnd(N,D);', ...
    'error(peel(lasterr));');

[num,den] = lmf2tf(N, D, tol);
num = mat2pol(num); den = mat2pol(den);
num.v = Var; den.v = Var;
num.h = h; den.h = h;

%end .. lmf2rat
