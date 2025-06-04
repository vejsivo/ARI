function C = plus(A,B,tol)
%PLUS    Add two-sided polynomials
%              C = A + B
%
% C = A + B or C = PLUS(A,B) adds the tsp matrices A and B
% with zeroing according to the global zeroing tolerance.
%
% C = PLUS(A,B,TOL) works with zeroing specified by the input
% relative tolerance TOL.
%
% See also TSP/MINUS, TSP/UPLUS.

%       Author: J.Jezek 11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $   $Date: 29-Sep-1999  13:00:00  $
%                          $Date: 22-May-2000  12:00:00  $
%                          $Date: 31-Oct-2000  13:00:00  $
%                          $Date; 24-Jan-2002            $

global PGLOBAL;

na = nargin;
if na==2,
   tol = PGLOBAL.ZEROING;
elseif na<2,
   error('Not enough input arguments.');
end;

eval('A = tsp(A);','error(''Invalid 1st argument.'');');
eval('B = tsp(B);', ...
   'eval(''B = mdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(B,'mdf'),
   C = 0;
   eval('C = plus(A,B,tol);','error(peel(lasterr));');
   eval('C = tsp(C);',';');
   return;
end;   

[th,Ch,A,B] = testht(A,B);
if th==0,
   warning('Inconsistent sampling periods.');
end;
oo = min(A.o,B.o);
AA = shift(A.p,A.o-oo,'z'); BB = shift(B.p,B.o-oo,'z');
CC = 0;
eval('CC = plus(AA,BB,tol);','error(peel(lasterr));');
C = tsp(CC); C.h = Ch;
C = shift(C,oo);

%end .. @tsp/plus
