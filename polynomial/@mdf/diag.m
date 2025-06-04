function D = diag(A,k,tol)
%DIAG  Diagonal matrix or diagonal of matrix-den fraction
%
% Let Av be a matrix-den vector fraction with N components (i.e.
% a matrix fraction 1-by-N or N-by-1). Then
%    DIAG(Av,K)  
% is a square matrix-den fraction of dimensions
% (N+ABS(K))-by-(N+ABS(K)) with the elements of Av on the K-th
% diagonal. K = 0 is the main diagonal, K > 0 is above the main
% diagonal, and K < 0 is below the main diagonal.
% DIAG(Av) is the same as DIAG(Av,0).
% 
% Let Am be a matrix-den fraction (other than a vector). Then
%    DIAG(Am,K)
% is a column vector matrix-den fraction formed from the elements
% of the K-th diagonal of Am. DIAG(Am) is the main diagonal of Am.
%
% DIAG(DIAG(Am)) is a diagonal matrix.
%
% See also MDF/TRIL, MDF/TRIU.

%       Author:  J. Jezek  04-Jan-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 30_Sep-2002 $
%                     $ Date 13-Oct-2002  double(logical)  $
%                     $ Date 14-Oct-2002 $
%                     $ Date 28-Feb-2003 $

global PGLOBAL;

if nargin==1, k=0;
elseif ~isa(k,'double'),
   error('Invalid 2nd argument.');
end;
if nargin==3,
   if ~isa(tol,'double'),
      error('Invalid 3rd argument.');
   end;
end;

NN = 0;
eval('NN = diag(A.frac.num,k);','error(peel(lasterr));');
DD = diag(A.frac.den,k) + double(~diag(ones(A.frac.s),k));
D = mdf(NN,DD);
if strcmp(A.frac.c,'cop'), props(D,'cop',A.frac.tc);
end;
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

%end .. @mdf/diag
