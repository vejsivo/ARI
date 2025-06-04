function D = diag(A,k,tol)
%DIAG  Diagonal matrix or diagonal of right-den fraction
%
% Let Av be a right-den vector fraction with N components. Then
%    DIAG(Av,K)  
% is a square right-den fraction of dimensions
% (N+ABS(K))-by-(N+ABS(K)) with the elements of Av on the K-th
% diagonal. K = 0 is the main diagonal, K > 0 is above the main
% diagonal, and K < 0 is below the main diagonal.
% DIAG(Av) is the same as DIAG(Av,0).
% 
% Let Am be a right-den matrix fraction. Then
%    DIAG(Am,K)
% is a column vector right-den fraction formed from the elements of
% the K-th diagonal of Am. DIAG(Am) is the main diagonal of Am.
%
% DIAG(DIAG(Am)) is a diagonal matrix.
%
% DIAG(A,K,TOL)  or  DIAG(A,[],tol)  works with zeroing
% tolerance specified by TOL.
%
% See also RDF/TRIL, RDF/TRIU.

%       Author:  J. Jezek  07-Feb-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 30-Sep-2002 $

global PGLOBAL;

if nargin==1 | isempty(k), k=0;
elseif ~isa(k,'double'),
   error('Invalid 2nd argument.');
end;
if nargin<=2,
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid 3rd argument.');
end;

eval('A = sdf(A,tol);','error(peel(lasterr));');
D = rdf(diag(A,k));

if strcmp(PGLOBAL.COPRIME,'cop'), D = coprime(D,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), D = reduce(D,tol);
else D = smreduce(D);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), D = defract(D);
end;

%end .. @rdf/diag
