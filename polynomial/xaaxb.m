function X = xaaxb(A,B,varargin)
%XAAXB   Symmetric equation solver
%
% The command 
%    X = XAAXB(A,B)
% solves the bilateral symmetric matrix equation
%    XA' + AX' = B
% where A is a square matrix
% and B is a symmetric (Hermitian) matrix
%    B = B' .
% The solution is not unique.
%
% This macro exists only for completeness. For more details
% and for possible further arguments, see POL/XAAXB.

%      Author: J. Jezek 06-Aug-2001
%      Copyroght(c) 2001 by Polyx, Ltd.

if nargin<2,
   error('Not enough input arguments.');
end;

eval('X = xaaxb(pol(A),pol(B),varargin{:});', ...
   'error(peel(lasterr));');

%end .. xaaxb
