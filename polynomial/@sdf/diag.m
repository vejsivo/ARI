function D = diag(A,k,tol)
%DIAG   Diagonal matrix or diagonal of scalar-den fraction
%
% Let Av be a scalar-den fraction 1-by-N or N-by-1. Then  DIAG(Av,K)  
% is a square scalar-den fraction, having dimensions
% (N+ABS(K))-by-(N+ABS(K)) with the elements of Av on the K-th
% diagonal. K = 0 is the main diagonal, K > 0 is above the main
% diagonal, and K < 0 is below the main diagonal.
% DIAG(Av) is the same as DIAG(Av,0).
% 
% Let Am be a scalar-den fraction other than 1-by-N or N-by-1. Then
%    DIAG(Am,K)
% is a scalar-den fraction, which is a column vector, formed from the
% elements of the K-th diagonal of Am. DIAG(Am) is the main
% diagonal of Am.
%
% DIAG(DIAG(Am)) is a diagonal matrix.
%
% See also SDF/TRIL, SDF/TRIU.

%       Author:  J. Jezek  28-Jan-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $
%                     $ Date 28-Feb-2003 $

global PGLOBAL;

if nargin==1, k=0;
elseif ~isa(k,'double'),
   error('Invalid 2nd argument.');
end;
if nargin==3,
   if~isa(tol,'double'),
      error('Invalid 3rd argument.');
   end;
end;

NN = 0;
eval('NN = diag(A.frac.num,k);','error(peel(lasterr));');
D = sdf(NN,A.frac.den);

if strcmp(A.frac.r,'red'), props(D,'red');
end;
if strcmp(A.frac.p,'prop'), props(D,'prop',A.frac.tp);
end;

if strcmp(PGLOBAL.COPRIME,'cop'), D = coprime(D);
end;
if strcmp(PGLOBAL.REDUCE,'red'), D = reduce(D);
else D = smreduce(D);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), D = defract(D);
end;

%end .. @sdf/diag
