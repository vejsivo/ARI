function [num,den] = rmf2rat(N,D,tol)
%RMF2RAT  Convert a right matrix fraction to rational form
%
% Given two polynomial matrices N and D, where D is square and
% nonsingular with the same number of columns as N, the function
%    [NUM,DEN] = RMF2RAT(N,D[,TOL])
% returns the transfer function H = N * D^{-1} in the form of two
% polynomial matrices NUM and DEN, containing the numerators
% and denominators, respectively, of the entries of the transfer 
% function matrix H.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: LMF2RAT, RAT2RMF.

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
eval('[N,D] = testdnd(N,D,''r'');', ...
   'error(peel(lasterr));');
Var = ''; h = 0;
eval('[Var,h,N,D] = testvhnd(N,D);', ...
   'error(peel(lasterr));');

[num,den] = rmf2tf(N, D, tol);
num = mat2pol(num); den = mat2pol(den);
num.v = Var; den.v = Var;
num.h = h; den.h = h;

%end .. rmf2rat
