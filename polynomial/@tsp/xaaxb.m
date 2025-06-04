function X = xaaxb(A,B,varargin)
%XAAXB  Symmetric polynomial equation solver
%
% The commmand
%    X = XAAXB(A,B) 
% solves the bilateral symmetric matrix polynomial equation
%    X(z^-1)A'(z) + A(z^-1)X'(z) = B(z)
% where A is a polynomial in z,
%       B is a para-Hermitian tsp matrix, that is,
%    B(z) = B'(z^-1)  .
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
% The commmand
%    XAAXB(A,B,'syl') 
% solves the equation with the Sylvester matrix method. This is the 
% default method. 
%
% The commmand
%    XAAXB(A,B,'tri') 
% returns a solution with a lower-triangular absolute coefficient matrix.  
%
% The command
%    XAAXB(A,B,'red') 
% solves the equation with polynomial reductions, a version of the Euclidean 
% division algorithm for polynomials. 
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TSP/AXXAB, TSP/AXYAB.

%        Author:  J. Jezek, August 11, 1999.
%        Copyright (c) by Polyx, Ltd.

ni = nargin;
if ni<2,
   error('Not enough input arguments.');
end;
eval('A = pol(A); B = tsp(B);','error(peel(lasterr));');

if ni>2,
   if isa(varargin{1},'cell'),
      varargin = varargin{1};
   end;
else
   varargin = [];
end;

lv = length(varargin);
if lv>=1,
   for i = 1:lv;
      arg = varargin{i};
      if ~isempty(arg),
         if ~isa(arg,'char') & ~isa(arg,'double'),
            error(['Invalid ',nth(i+1),' argument.']);
         end;
      end;
   end;
end;

A = conj(A.');
eval('X = axxab(A,B,varargin);','error(peel(lasterr));');
X = conj(X.');

%end .. @tsp/xaaxb

