function C = mtimes(A,B,tol)
%MTIMES (*)  Matrix multiply two-sided polynomials
%           C = A*B
%
% C = A*B or C = MTIMES(A,B) is the matrix product of the two-sided
% polynomial matrices A AND B. Any scalar tsp may multiply anything.
% Otherwise, the number of columns of A must equal the number of
% rows of B. Runs with zeroing using global zeroing tolerance. 
%
% C = MTIMES(A,B,TOL) works with zeroing specified by the input 
% relative tolerance TOL. 
%
% See also TSP/TIMES.

%       Author: J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 22-May-2000  12:00:00  $
%                         $Date: 31-Oct-2000  13:00:00  $
%                         $Date: 24-Jan-2002            $

global PGLOBAL;

na=nargin;
if na==2,
   tol=PGLOBAL.ZEROING;
elseif na==3,
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
else
   error('Not enough input arguments.');
end;

eval('A = tsp(A);','error(''Invalid 1st argument.'');');
eval('B = tsp(B);', ...
   'eval(''B = sdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(B,'sdf'),
   C = 0;
   eval('C = mtimes(A,B,tol);','error(peel(lasterr));');
   eval('C = tsp(C);',';');
   return;
end;

[th,Ch,A,B] = testht(A,B);
if th==0,
   warning('Inconsistent sampling periods.');
end;

CC = 0;
eval('CC = mtimes(A.p,B.p,tol);', 'error(peel(lasterr));');
C = tsp(CC); C.h = Ch;
C = shift(C, A.o + B.o);

%end .. @tsp/mtimes
