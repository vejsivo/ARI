function [X,Y,R,S] = axbyc(A,B,C,varargin)
%AXBYC  Matrix polynomial equation solver
%
% The command
%    [X0,Y0] = AXBYC(A,B) 
% finds a particular solution (X0,Y0) of the matrix polynomial Diophantine 
% equation        
%      AX + BY = C
% If no polynomial solution exists then all the entries of X0, Y0 are set 
% equal to NaN.
% 
% The command
%    [X0,Y0] = AXBYC(A,B,C,DEGREE) 
% seeks a solution [X0; Y0] of degree DEGREE. If DEGREE is not specified
% then a solution of minimum overall degree is computed by an iterative 
% scheme. If DEGREE is negative then the function directly computes an
% upper bound degree solution.
%
% The command
%    [X0,Y0] = AXBYC(A,B,C,'minx') 
% seeks a solution of minimum degree in X0. The commmand
%    [X0,Y0] = AXBYC(A,B,C,'miny') 
% seeks a solution of minimum degree in Y0. 
%
% The command
%    AXBYC(A,B,C,'syl') 
% solves the equation through the Sylvester matrix method. This is the 
% default method.
%
% The command
%    [X0,Y0,R,S] = AXBYC(A,B,C) 
% additionally computes the right null-space of [A B] so that all the 
% solutions to AX + BY = C may be parametrized as
%       X = X0 + RT
%       Y = Y0 + ST
% where T is an arbitrary polynomial matrix.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: POL/NULL, AXBY0, AXYBC.

%    Author: D. Henrion, September 14, 1998.
%    Modified by  J. Jezek, 30-May-2000.
%    Updated to 3.0 by D. Henrion, September 18, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.

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

if (rA ~= rB) | (rA ~= rC),
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
 D = vertcat(B.',A.');  
% D = [B.'; A.'];
else
 D = vertcat(A.',B.');
% D = [A.'; B.'];
end;

if nargout >= 3,
   Z = 0; L = 0;
   eval('[Z,L] = xab(D,C.'',varargin);', ...
      'error(peel(lasterr));');
else
   Z = 0;
   eval('Z = xab(D,C.'',varargin);', ...
      'error(peel(lasterr));');
end;

if any(any(isnan(Z))),
 X = pol(NaN*ones(cA, cC));
 Y = pol(NaN*ones(cB, cC));
elseif isempty(Z),
 X = pol(zeros(cA, 0));
 Y = pol(zeros(cB, 0));
elseif flip,
 Y = Z(:, 1:cB).';
 X = Z(:, cB+1:cA+cB).';
else
 X = Z(:, 1:cA).';
 Y = Z(:, cA+1:cA+cB).';
end;
X.v = Xv; Y.v = Xv;
X.h = Xh; Y.h = Xh; 

if nargout >= 3,
 if isempty(L),
  R = pol(zeros(cA, 0));
  S = pol(zeros(cB, 0));
 elseif flip,
  S = L(:, 1:cB).';
  R = L(:, cB+1:cA+cB).';
 else
  R = L(:, 1:cA).';
  S = L(:, cA+1:cA+cB).';
 end;
 R.v = Xv; S.v = Xv;
 R.h = Xh; S.h = Xh;
end;

%end .. axbyc
