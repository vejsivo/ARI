function [X,K] = axbc(A,B,C,varargin)
%AXBC  Matrix polynomial equation solver
%
% The command
%    X0 = AXBC(A,B,C) 
% finds a particular solution of the linear matrix polynomial equation        
%    AXB = C
% If no polynomial solution exists then all the entries in X0 are set 
% equal to NaN.
%
% The command
%    X0 = AXBC(A,B,C,DEGREE) 
% seeks a solution X0 of degree DEGREE. If DEGREE is not specified
% then a solution of minimum overall degree is computed by an iterative 
% scheme. If DEGREE is negative, then the function directly computes
% an upper bound degree solution.
%
% The command
%    X0 = AXBC(A,B,C,'sqz') 
% seeks a solution X0 with 'squeezed' entry degrees. If N is the nullity of 
% KRON(B, A.') then the N downmost entry degrees of the column vector X0(:) 
% are minimized, at the expense of increasing the degrees of the other entries.
%
% If DEGREES is an array of the same size as X0 such that DEGREES(i,j) = 1 and 
% DEGREES(k,l) = 0 then the function X0 = AXBC(A,B,C,'sqz',DEGREES)
% attempts to minimize the (i,j)-th entry degree in X0, provided that
% the (k,l)-th entry degree may increase.
%
% The command
%    AXBC(A,B,C,'syl') 
% solves the equation through the Sylvester matrix method. This is the default 
% method.
%
% The command
%    [X0, K] = AXBC(A,B,C) 
% additionally returns a cell array K = {K1, K2, ..}, where the Ki are
% solutions to the homogeneous equation AKB = 0. All the solutions to
% the equation AXB = C may be parametrized as
%       X = X0 + T1*K1 + T2*K2 + .. + Tp*Kp
% where the Ti are arbitrary scalar polynomials.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: POL/NULL, AXB, XAB.

%    Author: D. Henrion, May 27, 1998.
%    Modified by: J.Jezek, May 24, 2000.
%    Updated to 3.0 by D. Henrion, September 18, 2000.
%    Copyright 1998 by Polyx, Ltd.

if nargin < 3,
 error('Not enough input arguments.');
end;
eval('A = pol(A); B = pol(B); C = pol(C);', ...
   'error(peel(lasterr));');

[tv,Xv,A,B,C] = testvp3(A,B,C);
if tv==2,
   error('Inconsistent variables.');
elseif tv==0,
   warning('Inconsistent variables.');
end;

[th,Xh,A,B,C] = testhp3(A,B,C,Xv);
if th==0,
   warning('Inconsistent sampling periods.');
end;

if any(any(isnan(A))) | any(any(isnan(B))) | any(any(isnan(C))) | ...
      any(any(isinf(A))) | any(any(isinf(B))) | any(any(isinf(C))),
   error('Polynomial is not finite.');
end;

[rA, cA] = size(A);
[rB, cB] = size(B);
[rC, cC] = size(C); 

% Handle scalar A or B.

if (rA == 1) & (cA == 1), A = eye(rC) * A; rA = rC; cA = rC;
end;
if (rB == 1) & (cB == 1), B = eye(cC) * B; rB = cC; cB = cC;
end;

if (rA ~= rC) | (cB ~= cC),
 error('Matrices of incompatible dimensions.');
end;

% Rearrange DEGREES array.

argin = {};
lv = length(varargin);
if lv>0,
 for i = 1:lv,
  arg = varargin{i}; 
  if ~isempty(arg) & isa(arg, 'double') & any(size(arg) - [1 1]), % matrix
   if any(sort(size(arg)) - sort([cA rB])),
    error(['Invalid degree index array; must be ' int2str(cA) 'x' int2str(rB) '.']);
   else,
    arg = arg(:);
   end;
  end;
  argin{i} = arg;
 end;
end;

% Use Kronecker product to call XAB.

A = kron(B, A.');
B = C(:).';

if nargout < 2,
 eval('Y = xab(A,B,argin);', 'error(peel(lasterr));');
else
 eval('[Y,L] = xab(A,B,argin);', 'error(peel(lasterr));');
end;

if any(any(isnan(Y))),
 X = pol(NaN*ones(cA, rB));
elseif isempty(Y),
 X = pol(zeros(0, rB));
else
 X = pol(zeros(cA, rB));
 X(:) = Y.';
end;
X.v = Xv; X.h = Xh;

if nargout >= 2,
 K = {};
 Y = pol(zeros(cA, rB));
 if ~isempty(L),
  for i = 1:size(L, 1),
   Y(:) = L(i, :).'; Y.v = Xv; Y.h = Xh;
   K{i} = Y;
  end;
 end;
end;

%end .. axbc
