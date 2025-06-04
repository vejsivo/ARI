function [X,Y] = axyab(A,B,varargin)
%AXXAB   Non-symmetric polynomial equation solver
%
% The command
%    [X,Y] = AXYAB(A,B)
% solves the scalar polynomial equation
%    A'X + Y'A = B .
% A is a stable polynomial in variable 'z' or 'z^-1',
% B is a two-sided polynomial, and
%    MAX(ABS(DEG(B)),ABS(TDEG(B))) <= DEG(A) .
%
% Note: not to be misshandled with the other routines,
% the second argument B should be explicitly of class 'tsp'.
%
% The solution X is a polynomial in the same variable as A.
% If the is no solution of such a degree then X is set to NaN.
%
% The commmand
%    AXYAB(A,B,'syl') 
% uses the Sylvester matrix algorithm. This is the default method.
%
% The commmand
%    AXYAB(A,B,'red') 
% uses the polynomial reduction algorithm, a version of the Euclidean 
% division algorithm for polynomials.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: TSP/AXXAB.

%    Author:  J. Jezek, August 11, 1999
%    Copyright (c) 1999 by Polyx, Ltd.
%    $ Revision $   $ Date 24-May-2000 $

if nargin<2,
   error('Not enough input arguments.');
end;
eval('A = pol(A); B = tsp(B);','error(peel(lasterr));');

ni = length(varargin);
if ni>0,
   for i = 1:ni;
      arg = varargin{i};
      if ~isa(arg,'char') & ~isa(arg,'double'),
         error(['Invalid ',nth(i+2),' argument.']);
      end;
   end;
end;

var = A.v;
if isempty(var),
   var = 'z',
elseif strcmp(var,'z'),
elseif strcmp(var,'z^-1'),
   B = mirror(B);
else
   error('Invalid variable symbol in 1st argument; must be ''z'' or ''z^-1''.');
end;

if any([size(A),size(B)]~=1),
   error('Invalid arguments; must be scalar polynomials.');
end;

AA = tsp(A); [th,Xh,AA,B] = testht(AA,B); A = pol(AA);
if th==0,
   warning('Inconsistent sampling periods.');
end;

dA = deg(A); dB = deg(B); tdB = tdeg(B);
if ~isfinite(dA) | (isfinite(dB) & max(abs([dB,tdB])) > dA),
   error('Inconsistent degrees of polynomials.');
end;

Br = nneg(shift(B,dA)); Bl = pol(0);
symbol(A,'q'); symbol(Bl,'q'); symbol(Br,'q');
eval('[X,Y] = axyab(A,Bl,Br,0,varargin{1:ni});','error(peel(lasterr));');
symbol(X,var); symbol(Y,var);
X.h = Xh; Y.h = Xh;

%end .. @tsp/axyab


 

