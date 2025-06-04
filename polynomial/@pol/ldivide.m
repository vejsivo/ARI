function R = ldivide(D,N)
%LDIVIDE(.\)  Element-wise left divide polynomials
%   R = D.\N  or  R = LDIVIDE(D,N)
%
% The command returns a matrix whose (i,j)-th entry is D(i,j)\N(i,j).
%
% If N is a polynomial matrix and D is a matrix of numbers then the
% result is a polynomial matrix. If D is a polynomial matrix then
% the result is a matrix-denominator fraction; all entries of D must
% be nonzero.
%
% Matrices D,N must have the same dimensions, unless one of them is
% scalar. In such a case, this scalar is divided into (or by)
% every entry of the second matrix.
%
% The variable symbols of polynomials  D,N  should be the same. When 
% not, a warning is issued and the symbols are changed to the 
% standard one. However, if one symbol is 'z' and the other 'z^-1'
% then the symbols play a role, no warning being issued, the resulting
% symbol is taken from N.
%
% See also MDF/MDF, POL/MLDIVIDE, POL/RDIVIDE.

%     Author:  J. Jezek   27-Dec-1999
%     Copyright(c) 1999 by Polyx, Ltd.
%     $ Revision 3.0 $  $ Date 21-Apr-2000 $
%                       $ Date 02-Aug-2000 $
%                       $ Date 24-Jan-2002 $
%                       $ Date 30-Sep-2002 $

global PGLOBAL;

if nargin<2,
   error('Not enough input arguments.');
end;
eval('N = pol(N);',...
   'eval(''N = mdf(N);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(D,'double'),
   D = D.^-1;
   eval('R = D.*N;','error(peel(lasterr));');
else
   eval('D = pol(D);','error(''Invalid 1st argument.'');');
   if isa(N,'mdf'),
      eval('R = ldivide(D,N);','error(peel(lasterr));');
   else
      eval('R = mdf(N,D);', 'error(peel(lasterr));');
      if strcmp(PGLOBAL.COPRIME,'cop'), R = coprime(R);
      end;
      if strcmp(PGLOBAL.REDUCE,'red'), R = reduce(R);
      else R = smreduce(R);
      end;
      if strcmp(PGLOBAL.DEFRACT,'defr'), R = defract(R);
      end;
   end;
end;
  
%end .. @pol/ldivide
