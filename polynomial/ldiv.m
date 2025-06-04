function [Q,R] = ldiv(N,D,tol)
%LDIV  Left polynomial matrix division
% 
% The command
%    [Q,R] = LDIV(N,D) 
% computes the polynomial matrix quotient Q and the polynomial 
% matrix remainder R such that
%      N = D*Q + R
% and the degree of R is strictly less than the degree of D.
% Moreover, if D is row reduced then the i-th row degree of R
% is strictly less than the i-th row degree of D.
% If D is nonsingular and row-reduced then the rational matrix
% D^(-1)*R is strictly proper and the matrices Q and R are unique.
%
% There may be no solution if D is singular. In this case all the 
% entries in Q and R are set to NaN.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: RDIV, ROWRED.

%    Author: D. Henrion,  May 13, 1998.
%    Copyright 1998 by Polyx, Ltd
%    Modified by D. Henrion, August 3, 1998. Bug with TOL.
%                J. Jezek,   Feb 28,   2003, error check

% The function is dual to RDIV.

if nargin <= 1,
 error('Not enough input arguments.');
elseif nargin == 2,
 eval('[Q, R] = rdiv(N.'', D.'');', 'error(peel(lasterr));');
else
 eval('[Q, R] = rdiv(N.'', D.'', tol);', 'error(peel(lasterr));'); 
end;

Q = Q.'; R = R.';

%end .. ldiv
