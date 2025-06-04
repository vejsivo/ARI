function [N,D] = rat2lmf(num, den, tol)
%RAT2LMF  Convert a rational matrix to a left polynomial matrix fraction
%
% Given two polynomial arrays NUM and DEN of the same dimensions,
% containing the numerators and the denominators of a rational 
% transfer matrix, respectively, the function
%    [N,D] = RAT2LMF(NUM, DEN [,TOL])
% converts the transfer matrix to the left coprime polynomial 
% matrix fraction D^-1 * N.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TF2LMF, LMF2RAT.

%    Author: R.C.W. Strijbos, November 13, 1998.
%    Modified by J. Jezek, Aug 20, 2001,  arg checking
%    Copyright 1998 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

switch nargin
case 2
   tol = PGLOBAL.ZEROING;
case 3
   if ~isa(tol,'double') | length(tol) > 1 | ...
         ~isreal(tol) | tol<0 | tol>1
      error('Invalid tolerance.')
   end
otherwise
   error('Not enough input arguments.');
end

eval('num = pol(num); den = pol(den);', ...
   'error(peel(lasterr));');

[td,num,den] = testdp(num,den);
if ~td,
   error('Matrices of inconsistent dimensions.');
end;

Var = ''; h = 0;
eval('[Var,h,num,den] = testvhnd(num,den);', ...
   'error(peel(lasterr));');

num = pol2mat(num); den = pol2mat(den);
[N,D] = tf2lmf(num,den,tol);
N.v = Var; D.v = Var;
N.h = h; D.h = h;

%end .. rat2lmf

