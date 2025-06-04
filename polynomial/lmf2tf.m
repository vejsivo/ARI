function [num,den] = lmf2tf(N,D,tol)
%LMF2TF  Converts a left matrix fraction to a Control System Toolbox transfer function
%
% Given two polynomial matrices N and D, where D is square and 
% nonsingular with the same number of columns as N, the function
%    [NUM,DEN] = LMF2TF(N, D [,TOL])
% returns the transfer function H = D^{-1} * N in Control System 
% Toolbox format.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: RMF2TF, TF2LMF.

%    Authors: D. Henrion, M. Sebek, R.C.W. Strijbos, October 23, 1998.
%    Updated to 3.0 by D. Henrion, July 25, 2000.
%                   by J. Jezek, July 23, 2001, arg checking
%    Copyright 1998-2000 by Polyx, Ltd.

%    The macro is dual to RMF2TF.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

switch nargin
case 2
   tol = PGLOBAL.ZEROING;
case 3
   if ~isempty(tol),
      if ~isa(tol,'double'),
         error('Invalid tolerance.')
      end;
   end;
otherwise
   error('Not enough input arguments.');
end

eval('N = pol(N); D = pol(D);', ...
   'error(peel(lasterr));');

% Call RMF2TF.

eval('[num,den] = rmf2tf(N.'', D.'', tol);', ...
   'error(peel(lasterr));');
num = pol2mat(mat2pol(num).');
den = pol2mat(mat2pol(den).');

%end .. lmf2tf
