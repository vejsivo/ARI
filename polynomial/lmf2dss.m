function [a,b,c,d,e] = lmf2dss(N,D,tol)
%LMF2DSS  Converts a left matrix fraction to a descriptor state space system
%
% Given two polynomial matrices N and D such that D is column reduced
% the command
%    [a,b,c,d,e] = LMF2DSS(N, D [,TOL])
% returns the descriptor state space realization (a,b,c,d,e) of the
% system with transfer matrix H, that is,
%     H(x) = D^-1(x) N(x) = c (xe-a)^-1 b + d,
% where x is 's' (continuous-time) or 'z' (discrete-time).
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: LMF2SS.

%    Author: R.C.W. Strijbos, December 11, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision $  $ Date 22-Jul-2001  J.Jezek, arg checking  $

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

eval('[A,B,C,Dm] = lmf2ss(N,D,tol);', ...
   'error(peel(lasterr));');
[a,b,c,d,e] = ss2dss(A,B,C,Dm,tol);

%end .. lmf2dss
