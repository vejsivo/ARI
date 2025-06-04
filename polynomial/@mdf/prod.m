function Rs = prod(R,dim,tol);
%PROD   Element-wise product of matrix-den fraction
%
% For vectors,  PROD(X) is the product of the elements of X. 
% For matrices, PROD(X) is a row vector with the product over each column. 
%
% PROD(X,DIM) works along the dimension DIM.
%
% PROD(X,DIM,TOL) or PROD(X,[],TOL) works with zeroing specified by the
% input tolerance TOL.

%       Author:  J. Jezek 03-Jan-200
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 20-Sep-2002 $
%                     $ Date 28-Feb-2003 $

global PGLOBAL;

eval('R = mdf(R);','error(peel(lasterr));');

ni = nargin;
if ni==1,
  dim = []; tol = PGLOBAL.ZEROING;
elseif ni==2 | isempty(tol),
  tol = PGLOBAL.ZEROING;
else
  if ~isa(tol,'double'),
     error('Invalid tolerance.');
  end;
end;
if ~isempty(dim) & ~isa(dim,'double'),
   error('Invalid dimension.');
end;

Rs = 0;
eval('Rs = mdf(prod(R.frac.num,dim,tol),prod(R.frac.den,dim,tol));', ...
   'error(peel(lasterr));');

if strcmp(R.frac.p,'prop'),
   props(Rs,'prop');
end;
if strcmp(R.frac.r,'red'),
   props(Rs,'red');
end;
if strcmp(PGLOBAL.COPRIME,'cop'), Rs = coprime(Rs,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Rs = reduce(Rs);
else Rs = smreduce(Rs);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Rs = defract(Rs);
end;

%end .. @mdf/prod
