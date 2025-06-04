function c = eq(A,B,tol)
%EQ  (==)  Test if two-sided polynomials equal
%           A == B
%
% The commmand
%    C = (A==B) 
% performs elementwise comparison between two-sided polynomial matrices 
% A and B with a tolerance activated through the global variable
% PGLOBAL.ZEROING. A and B must have the same dimensions unless 
% one is a scalar two-sided polynomial. The scalar is compared with
% every entry of the other matrix.
%
% The commmands
%    C = EQ(A,B) 
% works alike. The commmand
%    C = EQ(A,B,TOL) 
% works with tolerance specified by the input tolerance TOL.
%  
% See also: TSP/NE

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 05-Apr-2000  $
%                         $Date: 22-May-2000  $
%                         $Date: 29-Jan-2002  $

global PGLOBAL;

na = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else   
   if ~isa(tol,'double')
      error('Invalid tolerance.');
   end;
end;

eval('A = tsp(A);','error(''Invalid 1st argument.'');');
eval('B = tsp(B);', ...
   'eval(''B = mdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(B,'mdf'),
   c = 0;
   eval('c = eq(A,B,tol);','error(peel(lasterr));');
   return;
end;   


[th,Ch,A,B] = testht(A,B);
if th==0,
   warning('Inconsistent sampling periods.');
end;

oo = min(A.o,B.o);
AA = shift(A.p,A.o-oo,'z'); BB = shift(B.p,B.o-oo,'z');
eval('c = eq(AA,BB,tol);','error(peel(lasterr));');

%end .. @tsp/eq

   