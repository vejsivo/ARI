function D = diag(A,k,tol)
%DIAG  Diagonal matrix or diagonal of two-sided polynomial
%
% Let Av be a two-sided polynomial vector with N components. Then
%    DIAG(Av,K)  
% is a square two-sided polynomial matrix of dimensions
% (N+ABS(K))-by-(N+ABS(K)) with the elements of Av on the K-th
% diagonal. K = 0 is the main diagonal, K > 0 is above the main
% diagonal, and K < 0 is below the main diagonal.
% DIAG(Av) is the same as DIAG(Av,0).
% 
% Let Am be a two-sided polynomial matrix. Then
%    DIAG(Am,K)
% is a column two-sided polynomial vector formed from the elements of the K-th
% diagonal of Am. DIAG(Am) is the main diagonal of Am.
%
% DIAG(DIAG(Am)) is a diagonal matrix.
%
% See also TSP/TRIL, TSP/TRIU.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 29-Sep-1999  13:00:00  $
%                         $Date: 22-May-2000  12:00:00  $
%                         $Date: 28-Nov-2000  12:00:00  $
%                         $Date: 28-Feb-2003  $

if nargin==1, k=0;
elseif ~isa(k,'double'),
   error('Invalid 2nd argument.');
end;
if nargin==3,
   if ~isa(tol,'double'),
      error('Invalid 3rd argument.');
   end;
end;

PP = 0;
eval('PP = diag(A.p,k);','error(peel(lasterr));');
D = tsp(PP); D.h = A.h;
D = shift(D,A.o); 

%end .. @tsp/diag
