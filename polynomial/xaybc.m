function [X,Y,R,S] = xaybc(A,B,C,varargin)
%XAYBC  Matrix polynomial equation solver
%
% The commmand
%    [X0,Y0] = XAYBC(A,B) 
% finds a particular solution of the matrix polynomial Diophantine equation   
%      XA + YB = C
% If no polynomial solution exists then all the entries in X0, Y0 are set 
% equal to NaN.
%
% The commmand
%    [X0,Y0] = XAYBC(A,B,C,DEGREE) 
% seeks a solution [X0 Y0] of degree DEGREE. If DEGREE is not specified then 
% a solution of minimum overall degree is computed by an iterative scheme.
% If DEGREE is negative then the function directly computes an upper bound 
% degree solution.
%
% The commmand
%    [X0,Y0] = XAYBC(A,B,C,'minx') 
% seeks a solution of minimum degree in X0. The commmand
%    [X0,Y0] = XAYBC(A,B,C,'miny') 
% seeks a solution of minimum degree in Y0.
%
% The commmand
%    XAYBC(A,B,C,'syl') 
% solves the equation through the Sylvester matrix method. This is the default
% method.
%
% The commmand
%    [X0,Y0,R,S] = XAYBC(A,B,C) 
% additionally computes the left null-space of [A; B] so that all the solutions 
% to XA + YB = C may be parametrized as
%       X = X0 + TR
%       Y = Y0 + TS
% where T is an arbitrary polynomial matrix.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: POL/NULL, XAB.

%    Author: D. Henrion, September 14, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    Modified by J. Jezek, August 3, 2001.

if nargin < 3,
 error('Not enough input arguments.');
elseif nargin < 4,
 varargin = [];
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


if any(any(isnan(A))) | any(any(isnan(B))) | any(any(isnan(C))),
 error('Polynomial is not finite.');
elseif any(any(isinf(A))) | any(any(isinf(B))) | any(any(isinf(C))),
 error('Polynomial is not finite.');
end;

[rA, cA] = size(A);
[rB, cB] = size(B);
[rC, cC] = size(C); 

if (cA ~= cB) | (cA ~= cC),
 error('Matrices of inconsistent dimensions.');
end;

flip = 0;
lv = length(varargin);
if lv>0,
 for i = 1:lv,
  arg = varargin{i};
  if isa(arg, 'char'),
   if strcmp(arg, 'minx'),
    varargin{i} = 'sqz'; % degree squeezing in X
    flip = 1;
   elseif strcmp(arg, 'miny'),
    varargin{i} = 'sqz';
   end;
  end;
 end;  
end;

% swap A and B if minimum degree in X is required
if flip,
   % D = [B; A];
   D = vartcat(B,A);
else
   %D = [A; B];
   D = vertcat(A,B);
end;

if nargout >= 3,
   Z = 0; L = 0;
   eval('[Z,L] = xab(D,C,varargin);', ...
      'error(peel(lasterr));');
else
   Z = 0;
   eval('Z = xab(D,C,varargin);', ...
      'error(peel(lasterr));');
end;

if any(any(isnan(Z))),
 X = pol(NaN*ones(rC, rA));
 Y = pol(NaN*ones(rC, rB));
elseif isempty(Z),
 X = pol(zeros(0, rA));
 Y = pol(zeros(0, rB));
elseif flip,
 Y = Z(:, 1:rB);
 X = Z(:, rB+1:rA+rB);
else
 X = Z(:, 1:rA);
 Y = Z(:, rA+1:rA+rB);
end;
X.v = Xv; Y.v = Xv;
X.h = Xh; Y.h = Xh;

if nargout >=3,
 if isempty(L),
  R = pol(ones(0, rA));
  S = pol(ones(0, rB));
 elseif flip,
  S = L(:, 1:rB);
  R = L(:, rB+1:rA+rB);
 else
  R = L(:, 1:rA);
  S = L(:, rA+1:rA+rB); 
 end;
 R.v = Xv; S.v = Xv;
 R.h = Xh; S.h = Xh;
end;

%end .. xaybc
