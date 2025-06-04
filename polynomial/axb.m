function [X,K] = axb(A,B,varargin)
%AXB  Matrix polynomial equation solver
%
% The command
%    X0 = AXB(A,B) 
% finds a particular solution X0 of the linear matrix polynomial 
% equation        
%    AX = B
% If no polynomial solution exists then all the entries in X0 are 
% set equal to NaN.
%
% The command
%    X0 = AXB(A,B,DEGREE) 
% seeks a solution X0 of degree DEGREE. If DEGREE is not specified
% then a solution of minimum overall degree is computed by an iterative 
% scheme. If DEGREE is negative then the function directly computes
% an upper bound degree solution.
%
% The command
%    X0 = AXB(A,B,'sqz') 
% seeks a solution X0 with 'squeezed' row degrees. If N is the nullity 
% of A then the N downmost row degrees are minimized, at the expense of 
% increasing the degrees in the other rows. 

% If DEGREES is a vector of zeros and ones such that DEGREES(i) = 1 and 
% DEGREES(j) = 0 then the function 
%    X0 = AXB(A,B,'sqz',DEGREES)
% attempts to minimize the degree of the i-th row in X0, provided that
% the degree of the j-th row may increase.
%
% The command
%    AXB(A,B,'syl') 
% solves the equation through the Sylvester matrix method. This is the 
% default method.
%
% The command
%    [X0,K] = AXB(A,B) 
% also computes the right null-space of A so that all the solutions to 
% AX = B may be parametrized as
%       X = X0 + KT
% where T is an arbitrary polynomial matrix.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: POL/NULL, XAB.

%    Author: D. Henrion, May 27, 1998.
%    Updated to 3.0 by D. Henrion, September 18, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

% Call XAB with tranposed arguments.

if nargin < 2,
 error('Not enough input arguments.');
elseif nargin < 3,
 varargin = [];
end;

if nargout < 2,
 eval('X = xab(A.'', B.'', varargin).'';','error(peel(lasterr));');
else
 eval('[X, K] = xab(A.'', B.'', varargin);','error(peel(lasterr));');
 X = X.'; K = K.';
end;

%end .. axb
