function [D,rk,U,Ui] = rowred(A, varargin)
%ROWRED   Row-reduced form of a polynomial matrix
%
% The command
%    [D,RK,U,UI] = ROWRED(A[,METHOD][,TOL]) 
% brings the polynomial matrix A into row reduced form D = U*A, 
% where U is an unimodular matrix and UI its inverse. The rows 
% of D are arranged according to decreasing degrees.
%
% If the input matrix A does not have full row rank then it is
% reduced to the form D = [D~' 0]' with D~ row reduced.
%
% The number of non-zero rows in D is returned as the scalar RK.
%
% The default method 'ref' ,is by repeated extraction of factors. A
% second method 'bas' is based on the calculation of a minimal basis.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: COLRED.

%    Authors: D. Henrion, R.C.W. Strijbos, September 18, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision 3.0 $  $ Date 30-May-2000  J.Jezek $

% The function is dual to COLRED

if nargin<1
   error('Not enough input arguments.');
end;
switch nargout
   case {0,1,2,3}
      eval('[D,rk,U] = colred(A.'',varargin{:});',...
         'error(peel(lasterr));');
      D = D.'; U = U.';
   case 4
      eval('[D,rk,U,Ui] = colred(A.'',varargin{:});',...
         'error(peel(lasterr));');
      D = D.'; U = U.'; Ui = Ui.';
end

%end .. rowred
