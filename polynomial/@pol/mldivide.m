function X = mldivide(A,B)
%MLDIVIDE (\)   Matrix left divide polynomials
%          X = A\B  or  X = MLDIVIDE(A,B)
%
% If A is a matrix of numbers and and B is a polynomial matrix then
% the command X = A\B returns a polynomial matrix. If A is a polynomial
% matrix then the command creates a left-denominator fraction with
% denominator A and numerator B.
%
% The dimensions of matrices A,B must be compatible (n-by-n, n-by-m)
% unless one of the matrices is scalar, which can be divided into
% everything and by everything.
%
% If both A and B are polynomials then their variable symbols should be
% the same. When not, a warning is issued and both symbols are changed
% to the standard one. However, if one symbol is 'z' and the other 'z^-1'
% then the symbols play a role, no warning, the resulting symbol is
% taken from B.
%
% See also: AXB, LDIV, POL/LDIVIDE, POL/MRDIVIDE, LDF/LDF.

%    Author:  J. Jezek  19-Nov-1999
%    Copyright 1999 by Polyx, Ltd.
%    $Revision $  $Date 21-Apr-2000   $
%                 $Date 25-Jan-2002   $
%                 $Date 30-Sep-2002   $

global PGLOBAL;

if nargin<2,
   error('Not enough input arguments.');
end;

eval('B = pol(B);',...
   'eval(''B = sdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(A,'double') & ndims(A)==2,
   eval('A = inv(A);', 'error(peelf(lasterr));');
   eval('X = A*B;', 'error(peel(lasterr));');
else
   eval('A = pol(A);','error(''Invalid 1st argument.'');');
   if isa(B,'sdf'),
      eval('X = mldivide(A,B);','error(peel(lasterr));');
   else
      eval('X = ldf(A,B);', 'error(peel(lasterr));');
      if strcmp(PGLOBAL.COPRIME,'cop'), X = coprime(X);
      end;
      if strcmp(PGLOBAL.REDUCE,'red'), X = reduce(X);
      else X = smreduce(X);
      end;
      if strcmp(PGLOBAL.DEFRACT,'defr'), X = defract(X);
      end;
   end;
end;

%end .. @pol/mldivide






