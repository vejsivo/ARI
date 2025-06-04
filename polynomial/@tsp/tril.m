function Au = tril(A,k,tol)
%TRIL    Lower triangular part of two-sided polynomial
%
% TRIL(P) is the lower triangular part of X.
% TRIL(P,K) is the elements on and below the K-th diagonal
% of P .  K = 0 is the main diagonal, K > 0 is above the
% main diagonal and K < 0 is below the main diagonal.
%
% TRIL(P,K,TOL) works with zeroing specified by the input 
% tolerance TOL.
%
% See also TSP/TRIU.

%        Author:  J. Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.
%        $Revision: 3.0 $  $Date: 29-Sep-1999 $
%                          $Date: 24-May-2000 $
%                          $Date: 31-Oct-2000 $

global PGLOBAL;

ni = nargin;
if ni<3, tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double'),
      error('Invalid 3rd argument.');
   end;
end;
if ni<2, k = 0;
else
   if ~isa(k,'double'),
      error('Invalid 2nd argument.');
   end;
end;

PP = 0;
eval('PP = tril(A.p,k,tol);','error(peel(lasterr));');
Au = tsp(PP); Au.h = A.h;
Au = shift(Au,A.o);

%end .. @tsp/tril
