function A = sum(R,dim,tol);
%SUM   Element-wise sum of scalar-den fraction
% 
% For a scalar-den fraction R which is a vector (i.e., a matrix
% having one row or one column), SUM(R) is the sum of the
% elements of R. 
%
% For other scalar-den fractions, SUM(R) is a row vector with
% the sums over each column. 
%
% SUM(R,DIM) sums along the dimension DIM, where DIM is 1 or 2.
%
% SUM(R,DIM,TOL) or SUM(R,[],TOL) works with zeroing specified by 
% the input tolerance TOL.

%      Author:  J. Jezek  28-Jan-2000
%      Copyright (c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 26-Apr-2000 $
%                    $ Date 06-Nov-2000 $
%                    $ Date 30-Sep-2002 $
%                    $ Date 14-Oct-2002 $
%                    $ Date 28-Feb-2003 $

global PGLOBAL;

ni = nargin;
if ni==1,
  dim = []; tol = PGLOBAL.ZEROING;
elseif ni==2 | isempty(tol),
  tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;
if ~isempty(dim) & (~isa(dim,'double') | (dim~=1 & dim~=2)),
   error('Invalid dimension.');
end;

NN = 0;
eval('NN = sum(R.frac.num,dim,tol);','error(peel(lasterr));');
A = sdf(NN,R.frac.den);

if strcmp(R.frac.p,'prop'),
   props(A,'prop');
end;
if strcmp(PGLOBAL.COPRIME,'cop'),
   A = coprime(A,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   A = reduce(A);
else
   A = smreduce(A);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   A = defract(A);
end;

%end .. @sdf/sum
