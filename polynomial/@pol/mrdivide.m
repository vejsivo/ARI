function X = mrdivide(B,A)
%MRDIVIDE (/)   Matrix right divide polynomials
%          X = B/A  or  X = MRDIVIDE(B,A)
%
% If B is a polynomial matrix and A is a matrix of numbers then
% the command  X = B/A  returns a polynomial matrix. If A is polynomial
% matrix then the command creates a right-denominator fraction with
% numerator B and denominator A.
%
% The dimensions of matrices  B,A  must be compatible  (m-by-n,
% n-by-n)  unless one of them is scalar, which can be divided
% by everything and into everything.
%
% If both B and A are polynomials then their variable symbols should
% be the same. When not, a warning is issued and both symbols are
% changed to the standard one. However, if one symbol is 'z' and the
% other 'z^-1' then the symbols play a role, no warning,
% the resulting symbol taken from B.
%
% See also: XAB, RDIV, POL/RDIVIDE, POL/MLDIVIDE, RDF/RDF.

%    Author:  J. Jezek  18-Nov-1999
%    Copyright 1999 by Polyx, Ltd.
%    $Revision$  $Date  21-Apr-2000  $
%                $Date  02-Aug-2000  $
%                $Date  25-Jan-2002  $
%                $Date  30-Sep-2002  $

global PGLOBAL;

if nargin<2,
   error('Not enough input arguments.');
end;

eval('B = pol(B);', 'error(''Invalid 1st argument.'');');
if isa(A,'double') & ndims(A)==2,
   eval('A = inv(A);', 'error(peelf(lasterr));');
   eval('X = B*A;', 'error(peel(lasterr));');
else
   eval('A = pol(A);',...
      'eval(''A = sdf(A);'',''error(''''Invalid 2nd argument.'''');'');');
   if isa(A,'sdf'),
      eval('X = mrdivide(B,A);','error(peel(lasterr));');
   else
      eval('X = rdf(B,A);', 'error(peel(lasterr));');
      if strcmp(PGLOBAL.COPRIME,'cop'), X = coprime(X);
      end;
      if strcmp(PGLOBAL.REDUCE,'red'), X = reduce(X);
      else X = smreduce(X);
      end;
      if strcmp(PGLOBAL.DEFRACT,'defr'), X = defract(X);
      end;
   end;
end;

%end .. @pol/mrdivide

