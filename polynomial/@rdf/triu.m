function Au = triu(A,k,tol)
%TRIU    Upper triangular part of right-den fraction
%
% TRIU(A) is the upper triangular part of X.
% TRIU(A,K) is the elements on and above the K-th diagonal
% of A.  K = 0 is the main diagonal, K > 0 is above the
% main diagonal and K < 0 is below the main diagonal.
%
% TRIU(A,K,TOL) works with zeroing specified by the input 
% tolerance TOL.
%
% See also RDF/TRIL.

%        Author:  J. Jezek  07-Feb-2000
%        Copyright (c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 25-Apr-2000 $
%                      $ Date 30-Sep-2002 $

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

eval('A = sdf(A,tol);','error(peel(lasterr));');
Au = rdf(triu(A,k));

if strcmp(PGLOBAL.COPRIME,'cop'), Au = coprime(Au);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Au = reduce(Au);
else Au = smreduce(Au);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Au = defract(Au);
end;

%end .. @rdf/triu
