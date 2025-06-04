function [D,rk,U,Ui] = colred(A,varargin);
%COLRED   Column-reduced form of a two-sided polynomial matrix
%
% The command
%    [D,RK,U,UI] = COLRED(A[,METHOD][,TOL])
% brings the two-sided polynomial matrix A into column reduced
% form  D = A*U  where U is a unimodular polynomial matrix and
% UI its inverse. If the input matrix A does not have full
% column rank then it is reduce to the form  D = [D~ 0] with
% column reduced. The number of nonzero column in D is returned
% as the scalar RK.
%
% If the input matrix A is square then the operation is modified:
% D is such that  D = U'*A*U . If A is para-Hermitian then D
% is also.  If A does not have full column rank then it is reduced
% to the form   [D~ 0
%                0  0]  with RK nonzero rows and columns.
%
% The additional input arguments METHOD and TOL are as for COLRED.
%
% See also COLRED.

%      Author:  J. Jezek  11-Mar-2003
%      Copyright(c) 2003 by Polyx, Ltd.

if ~isa(A,'tsp'),
   error('Some argument but not 1st is invalidly tsp.');
end;

PP = 0; rk = 0; U = 0; Ui = 0;
if nargout==4,
   eval('[PP,rk,U,Ui] = colred(A.p,varargin{:});', ...
      'error(peel(lasterr));');
else
   eval('[PP,rk,U] = colred(A.p,varargin{:});', ...
      'error(peel(lasterr));');
end;
D = tsp(PP); D.h = A.h;
D = shift(D,A.o);
[As1,As2] = size(A);
if As1==As2,
   D = U'*D;
end;

%end .. @tsp/colred

   