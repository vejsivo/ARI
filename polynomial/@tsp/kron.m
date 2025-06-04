function C = kron(A,B,tol)
%KRON     Kronecker tensor product of two-sided polynomials
%
% KRON(A,B) is the Kronecker tensor product of A and B.
% The result is a large matrix formed by taking all possible
% products between the elements of A and those of B.  
%
% KRON(A,B,TOL) works with zeroing specified by the input 
% relative tolerance TOL.
%
% See also TSP/TIMES, TSP/MTIMES.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 22-May-2000  12:00:00  $
%                         $Date: 31-Oct-2000  13:00:00  $
%                         $Date: 25-Jan-2002            $

global PGLOBAL;

na=nargin;
if na==2,
   tol=PGLOBAL.ZEROING;
elseif na==3,
   if ~isa(tol,'double'),
      error('Invalid 3rd argument.');
   end;
else
   error('Not enough input arguments.');
end;

eval('A = tsp(A);','error(''Invalid 1st argument.'');');
eval('B = tsp(B);', ...
   'eval(''B = mdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(B,'mdf'),
   C = 0;
   eval('C = kron(A,B,tol);','error(peel(lasterr));');
   eval('C = tsp(C);',';');
   return;
end;

[th,Ch,A,B] = testht(A,B);
if th==0,
   warning('Inconsistent sampling periods.');
end;

CC = 0;
eval('CC = kron(A.p,B.p,tol);', 'error(peel(lasterr));');
C = tsp(CC); C.h = Ch;
C = shift(C, A.o + B.o);

%end .. @tsp/kron
