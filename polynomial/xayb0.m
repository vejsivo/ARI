function [X,Y] = xayb0(A,B,tol);
%XAYB0    Solver of polynomial equation  X*A + Y*B = 0
%
% For polynomial matrices A and B,
% number of columns of A == number of columns of B,
% matrix A square and nonsingular,
% the command
%     [X,Y] = XAYB0(A,B)
% returns polynomial matrices X and Y satisfying the equation
%     X*A + Y*B = 0 ,
% number of rows of X == number of rows of Y,
% matrix Y square and nonsingular.
% 
% A tolerance may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also AXBY0.

%      Author:  J. Jezek  19-Nov-1999
%      Copyright(c) 1999 by Polyx Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni<2,
   error('Not enough input arguments.');
elseif ni==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double')
      error('Invalid tolerance.');
   end;
end;

eval('A = pol(A); B = pol(B);', 'error(peel(lasterr));');

eval('[X,Y] = axby0(A.'',B.'',tol);', 'error(peel(lasterr));');
X = X.'; Y = Y.';

%end .. xayb0
