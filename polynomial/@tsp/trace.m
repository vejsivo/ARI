function At = trace(A,tol)
%TRACE   Sum of the diagonal elements of two-sided polynomial
%   
% TRACE(P) is the sum of the diagonal elements of P
% with zeroing activated through the global variable
% PGLOBAL.ZEROING.
%
% TRACE(P,TOL) works with zeroing specified by the input 
% tolerance TOL.

%       Author: J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 24-May-2000  12:00:00  $
%                         $Date: 31-Oct-2000  13:00:00  $

global PGLOBAL;

if nargin<2, tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

PP = 0;
eval('PP = trace(A.p,tol);','error(peel(lasterr));');
At = tsp(PP); At.h = A.h;
At = shift(At,A.o);

%end .. @tsp/trace
