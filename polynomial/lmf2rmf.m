function [P,Q] = lmf2rmf(N,D,tol)
%LMF2RMF  Left-to-right matrix fraction conversion
%
% Given a nonsingular polynomial matrix D and a polynomial matrix N with 
% the same number of rows, the command
%    [P,Q] = LMF2RMF(N,D)
% computes a nonsingular polynomial matrix Q and a polynomial matrix P
% such that
%    D^{-1}*N = P*Q^{-1}
% Moreover, P and Q are right coprime and if D^{-1}*N is proper
% then Q is column reduced.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: RMF2LMF.

%    Authors: D. Henrion, H. Kwakernaak, S. Pejchova,  May 12, 1998.
%    Version 2.0 last modified by D. Henrion, August 25, 1999.
%    Updated to version 3.0 by D. Henrion, July 25, 2000.
%    Modified  by J. Jezek July 23, 2001, arg check, sampling periods
%    Copyright 1998-2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin < 2
   error('Not enough input arguments.');
end;

eval('N = pol(N);', 'error(peel(lasterr));');
eval('D = pol(D);', 'error(peel(lasterr));');
eval('[N,D] = testdnd(N,D,''l'');', 'error(peel(lasterr));');

[rN, cN] = size(N); [rD, cD] = size(D);
% Empty matrices: section added by D. Henrion, August 1, 2000
if (rD == 0) | (cN == 0), % empty matrices
   P = pol(zeros(rD,cN)); Q = pol(eye(cN));
   return;
end;

if nargin < 3,
   tol = [];
elseif ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

[tv,Var,N,D] = testvp(N,D);
if tv==2,
   [th,h,N,D] = testhp(N,D,Var);
   if ~th,
      warning('Inconsistent sampling periods.');
   end;
   if strcmp(D.v,'z^-1'),
      dg = deg(D);
      N = shift(N,dg); D = rev(D,dg);
      D.v = 'z';
   else
      dg = deg(N);
      N = rev(N,dg); D = shift(D,dg);
      N.v = 'z';
   end;
   Var = 'z';
elseif ~tv,
   warning('Inconsistent variables.');
end;

[th,h,N,D] = testhp(N,D,Var);
if ~th,
   warning('Inconsistent sampling periods.');
end;

me = 1; % relative tolerance for zeroing
if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else
 tolzero = tol * me;
end;

if issingular(D, tolzero),
  error('Denominator matrix is singular.');
end;

% Compute the right nullspace of [D -N]

PQ = null(horzcat(D,-N), tol);
PQ.h = h;

% Identify P and Q

P = pzer(PQ(1:cD, :), tolzero);
Q = pzer(PQ(cD+1:cD+cN, :), tolzero);

%end .. lmf2rmf
