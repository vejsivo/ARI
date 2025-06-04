function X = axxab(A,B,varargin)
%AXXAB  Symmetric polynomial equation solver
%
% The command
%    X = AXXAB(A,B)
% solves the bilateral symmetric matrix polynomial equation
%    A'(z^-1)X(z) + X'(z^-1)A(z) = B(z)
% where A is a square polynomial matrix in variable 'z' or 'z^-1',
% B is a para-Hermitian two-sided polynomial matrix, that is,
%    B(z) = B'(z^-1) .
%
% Note: not to be misshandled with the other routines,
% the second argument B should be explicitly of class 'tsp'.
%
% The solution X is a square polynomial matrix in the same
% variable as A, of degree DEG(X) = MAX(DEG(A), DEG(B)).
% If there is no solution of such a degree then all the entries 
% in X are set equal to NaN.
%
% The macro may be used with several modifiers:
%
% The command
%    AXXAB(A,B,'syl') 
% solves the equation with the Sylvester matrix method. This is the 
% default method.
%
% The command
%    AXXAB(A,B,'tri') 
% returns a solution with upper-triangular absolute coefficient matrix.
%
% The command
%    AXXAB(A,B,'red') 
% solves the equation with polynomial reductions, a version of the Euclidean 
% division algorithm for polynomials. 
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TSP/AXYAB.

%      Author:  J. Jezek, August 11, 1999.
%      Copyright 1999 by Polyx, Ltd.
%      $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                        $Date: 24-May-2000  12:00:00  $
%                        $Date: 06-Aug-2001  A(z),A(z^-1) $

ni = nargin;
if ni<2,
   error('Not enough input arguments.');
end;
eval('A = pol(A); B = tsp(B);','error(peel(lasterr));');

if isempty(A.v) | strcmp(A.v,'z'),
   var = 'z';
elseif strcmp(A.v,'z^-1'),
   var = 'z^-1';
   B = mirror(B);
else
   error('Invalid variable symbol in 1st argument; must be ''z'' or ''z^-1''.');
end;
if B~=B',
   error('Matrix is not para-Hermitian.');
end;

AA = tsp(A); [th,Xh,AA,B] = testht(AA,B); A = pol(AA);
if th==0,
   warning('Inconsistent sampling periods.');
end;

B = nneg(shift(B,B.r/2));
symbol(A,'q'); symbol(B,'q');

if ni>2,
   if isa(varargin{1},'cell'),
      varargin = varargin{1};
   end;
else
   varargin = [];
end;

lv = length(varargin);
if lv>0,
   for i = 1:lv;
      arg = varargin{i};
      if ~isempty(arg),
         if ~isa(arg,'char') & ~isa(arg,'double'),
            error(['Invalid ',nth(i+2),' argument.']);
         end;
      end;
   end;
end;

eval('X = axxab(A,B,varargin);','error(peel(lasterr));');
symbol(X,var); X.h = Xh;

%end .. @tsp/axxab
