function Au = triu(A,k,tol)
%TRIU    Upper triangular part of scalar-den fraction
%
% TRIU(A) is the upper triangular part of X.
% TRIU(A,K) is the elements on and above the K-th diagonal
% of A.  K = 0 is the main diagonal, K > 0 is above the
% main diagonal and K < 0 is below the main diagonal.
%
% TRIU(A,K,TOL) works with zeroing specified by the input 
% tolerance TOL.
%
% See also SDF/TRIL.

%        Author:  J. Jezek  28-Jan-2000
%        Copyright (c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 26-Apr-2000 $
%                      $ Date 06-Nov-2000 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 14-Oct-2002 $

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

NN = 0;
eval('NN = triu(A.frac.num,k);','error(peel(lasterr));');
Au = sdf(NN,A.frac.den);

if strcmp(A.frac.r,'red'), props(Au,'red');
end;
if strcmp(A.frac.p,'prop'), props(Au,'prop',A.frac.tp);
end;

if strcmp(PGLOBAL.COPRIME,'cop'), Au = coprime(Au);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Au = reduce(Au);
else Au = smreduce(Au);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Au = defract(Au);
end;

%end .. @sdf/triu
