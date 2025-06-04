function [N,D] = rat2rmf(num, den, tol)
%RAT2RMF  Convert a rational matrix to a right polynomial matrix fraction
%
% Given two polynomial arrays NUM and DEN of the same dimensions,
% containing the numerators and the denominators of the entries of
% a rational transfer matrix, respectively, the function
%    [N,D] = RAT2RMF(NUM,DEN[,TOL])
% converts the transfer matrix to the right coprime polynomial 
% matrix fraction  N * D^-1.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TF2RMF, RMF2RAT.

%    Author: R.C.W. Strijbos, November 13, 1998.
%    Modified by J. Jezek, Aug 20, 2001,  arg checkin
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
[N,D] = tf2rmf(num,den,tol);
N.v = Var; D.v = Var;
N.h = h; D.h = h;

%end .. rat2rmf

