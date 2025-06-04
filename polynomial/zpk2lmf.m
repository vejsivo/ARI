function [N,D] = zpk2lmf(Z,P,K,var,tol)
%ZPK2LMF  Convert from zero-pole-gain format to a left matrix fraction
%
% Given a rational transfer function H in zero-pole-gain format, that is,
% given cell-arrays Z, P such that Z{i,j} and P{i,j} are vectors of zeros
% and poles of H(i,j) and given the corresponding array of gains K,
% the function
%    [N,D] = ZPK2LMF(Z,P,K)
% returns a left coprime matrix fraction description of H = D^{-1}*N.
%
% A variable symbol may be specified as an additional input argument.
% Its default value is the global variable symbol.
% A tolerance parameter may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%    Author: D. Henrion, October 22, 1998.
%    Modified by D. Henrion, May 9, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

if nargin < 3,
 error('Not enough input arguments.');
elseif nargin < 4,
 var = PGLOBAL.VARIABLE;
 tol = PGLOBAL.ZEROING;
elseif nargin < 5,
 if ~isa(var,'char') & isa(var,'double'),
  tol = var;
  var = PGLOBAL.VARIABLE;
 else
  tol = PGLOBAL.ZEROING;
 end;
else
 if isa(tol,'char') & isa(var,'double'),
  tmp = var; var = tol; tol = tmp;
 end;
end;

if ~isa(tol,'double') | length(tol) > 1 | ...
      ~isreal(tol) | tol<0 | tol>1,
 error('Invalid tolerance.');
end;

if ~isa(var,'char'),
 error('Invalid variable symbol.');
end;

eval('N = root2pol(Z,K,tol);', 'error(peel(lasterr));');
eval('D = root2pol(P,[],tol);', 'error(peel(lasterr));');
N = pol2mat(N);
D = pol2mat(D);

[N,D] = tf2lmf(N, D, tol);
eval('N.v = var; D.v = var;', 'error(peel(lasterr));');

%end .. zpk2lmf
