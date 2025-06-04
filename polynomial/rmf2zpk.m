function [Z,P,K] = rmf2zpk(N,D,tol)
%RMF2ZPK  Convert a right matrix fraction to zero-pole-gain format
%
% Given two polynomial matrices N and D, where D is square and
% nonsingular with the same number of columns as N, the function
%     [Z,P,K] = RMF2ZPK(N, D)
% returns the transfer function H = N * D^{-1} in zero-pole-gain
% format. That is, the output arguments Z, P are cell-arrays
% such that Z{i,j} and P{i,j} are vectors containing the zeros and 
% poles of H(i,j), and the output argument K is the corresponding 
% array of gains.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%    Author: D. Henrion, October 22, 1998.
%    Updated to 3.0 by D. Henrion, July 25, 2000.
%                   by J. Jezek, July 23, 2001, arg check
%    Copyright 1998-2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

switch nargin
case 0,1
   error('Not enough input arguments.');
case 2
   tol = PGLOBAL.ZEROING;
case 3
   if ~isa(tol,'double'),
      error('Invalid tolerance.')
   end
end

eval('N = pol(N); D = pol(D);', ...
   'error(peel(lasterr));');

eval('[N,D] = rmf2tf(N,D,tol);', ...
   'error(peel(lasterr));');

[Z,KZ] = pol2root(mat2pol(N));
[P,KP] = pol2root(mat2pol(D));

K = KZ ./ KP;

%end .. rmf2zpk

