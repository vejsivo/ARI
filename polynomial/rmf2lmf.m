function [P,Q] = rmf2lmf(N,D,tol)
%RMF2LMF  Right-to-left matrix fraction conversion
%
% Given a non-singular polynomial matrix D and a polynomial matrix N 
% with the same number of columns, the command
%    [P,Q] = RMF2LMF(N,D)
% computes a non-singular polynomial matrix Q and a polynomial matrix P
% such that
%    N*D^{-1} = Q^{-1}*P
% Moreover, P and Q are left coprime and Q is row reduced.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: LMF2RMF.

%    Authors: D. Henrion, H. Kwakernaak, S. Pejchova,  May 12, 1998.
%    Version 2.0 last modified by D. Henrion, August 28, 1998.
%    Updated to version 3.0 by D. Henrion, July 25, 2000.
%    Last modified by D. Henrion, August 31, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

% The function is dual to LMF2RMF.

if nargin < 2
 error('Not enough input arguments.');
elseif nargin == 2,
 eval('[P,Q] = lmf2rmf(N.'', D.''); P = P.''; Q = Q.'';', ...
      'error(peel(lasterr));');
else
 eval('[P,Q] = lmf2rmf(N.'', D.'', tol); P = P.''; Q = Q.'';', ...
      'error(peel(lasterr));');
end

%end .. rmf2lmf

