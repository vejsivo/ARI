function [X,Y,R,S] = axybc(A,B,C,varargin)
%AXYBC  Matrix polynomial equation solver
%
% The command
%    [X0,Y0] = AXYBC(A,B,C) 
% finds a particular solution of the linear matrix polynomial equation    
%    AX + YB = C
% If no polynomial solution exists then all the entries in X0 are set 
% equal to NaN.
%
% The command
%    [X0,Y0] = AXYBC(A,B,C,DEGREE) 
% seeks a solution of degree DEGREE. If DEGREE is not specified then a 
% solution of minimum overall degree is computed by an iterative scheme.
% If DEGREE is negative then the function directly computes an upper bound 
% degree solution.
%
% The command
%    [X0,Y0] = AXYBC(A,B,C,'min') 
% seeks a solution such that the N downmost entry degrees of the column 
% vector Z = [X0(:); Y0(:)] are minimized. N is the nullity of the matrix 
% [KRON(I, A) KRON(B.', I)]. 

% If DEGREES is a vector of the same length as Z such that DEGREES(i) = 1 
% and DEGREES(j) = 0 then the function AXYBC(A,B,C,'min',DEGREES) attempts 
% to minimize the i-th entry degree in Z, provided that the k-th entry degree 
% may increase.
%
% The command
%    AXYBC(A,B,C,'syl') 
% solves the equation through the Sylvester matrix method. This is the default 
% method.
%
% The command
%    [X0,Y0,R,S] = AXYBC(A,B,C) 
% additionally returns two cell arrays R = {R1, R2, ..} and S = {S1, S2, ..} 
% such that the Ri and Si are solutions to the homogeneous equation
% AR + SB = 0. All the solutions to equation AX + YB = C may be
% parametrized as
%       X = X0 + T1*R1 + T2*R2 + .. + Tp*Rp
%       Y = Y0 + T1*S1 + T2*S2 + .. + Tp*Sp
% where the Ti are arbitrary scalar polynomials.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: POL/NULL, AXBYC.

%    Author: D. Henrion, June 3, 1998.
%    Modified by: J.Jezek, May 24, 2000.
%    Updated to 3.0 by D. Henrion, September 18, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

if nargin < 3,
 error('Not enough input arguments.');
end;

eval('A = pol(A); B = pol(B); C = pol(C);', ...
   'error(peel(lasterr));');

[tv,Xv,A,B,C] = testvp3(A,B,C);
if tv==0,
   warning('Inconsistent variables.');
end;

[th,Xh,A,B,C] = testhp3(A,B,C,Xv);
if th==0,
   warning('Inconsistent sampling periods.');
end;

if any(any(isnan(A))) | any(any(isnan(B))) | any(any(isnan(C))),
 error('Polynomial is not finite.');
elseif any(any(isinf(A))) | any(any(isinf(B))) | any(any(isinf(C))),
 error('Polynomial is not finite.');
end;

[rA, cA] = size(A);
[rB, cB] = size(B);
[rC, cC] = size(C); 

if (rA ~= rC) | (cB ~= cC),
 error('Matrices of inconsistent dimensions.');
end;

% Check DEGREES array.

ldeg = cA*cC+rC*rB;
argin = {};
for i = 1:length(varargin),
 arg = varargin{i};
 if ~isempty(arg) & isa(arg, 'double') & any(size(arg) - [1 1]), % matrix
  if any(sort(size(arg)) - [1 ldeg]),
   error(['The degree index vector must have ' int2str(ldeg) ' components.']);
  else,
   arg = arg(:);
  end;
 end;
 argin{i} = arg;
end;

% Use Kronecker product to call XAB.

A = [kron(eye(cB), A.'); kron(B, eye(rA))];
B = C(:).';
X = pol(zeros(cA, cC));
Y = pol(zeros(rC, rB));

if nargout < 3,

  Z = xab(A, B, varargin);

else
 
  [Z, T] = xab(A, B, varargin);

end;

if any(any(isnan(Z))),
 X = pol(NaN*ones(cA, cC));
 Y = pol(NaN*ones(rC, rB));
elseif isempty(Y),
 X = pol(zeros(cA, 0));
 Y = pol(zeros(0, rB));
else
 X = pol(zeros(cA, cC));
 Y = pol(zeros(rC, rB));
 X(:) = Z(1:cA*cC).';
 Y(:) = Z(cA*cC+1:cA*cC+rC*rB).';
end;
X.v = Xv; Y.v = Xv;
X.h = Xh; Y.h = Xh;

if nargout >=3,
 M = pol(zeros(cA, cC));
 N = pol(zeros(rC, rB));
 R = {}; S = {};
 if ~isempty(T),
  for i = 1:size(T, 1),
     M(:) = T(1:cA*cC).'; N(:) = T(cA*cC+1:cA*cC+rC*rB).';
     M.v = Xv; N.v = Xv; M.h = Xh; N.h = Xh;
     R{i} = M; S{i} = N;
  end;
 end;
end;

%end .. axybc
